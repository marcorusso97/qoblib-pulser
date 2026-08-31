#!/usr/bin/env python3
"""Eight-atom Pulser proof of concept for QOBLIB Steiner repairs.

Each benchmark case keeps all four candidates of the first two nets of a saved
``k=3, m=4`` LNS repair.  The resulting ``k=2, m=4`` problem has eight atoms,
16 one-hot selections, and a 256-dimensional Hilbert space.

The public Pulser ``AnalogDevice`` profile is not claimed to be a RUBY digital
twin.  The benchmark tests whether global Rydberg control can preserve the exact
repair optimum after a deterministic two-dimensional geometric fit.
"""
import argparse
import hashlib
import importlib.metadata
import itertools
import json
import math
import platform
import statistics
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pulser import Pulse, Register, Sequence
from pulser.backend.config import EmulationConfig
from pulser.backend.default_observables import StateResult
from pulser.devices import AnalogDevice
from pulser.waveforms import InterpolatedWaveform
from pulser_simulation import QutipBackendV2
from scipy.optimize import least_squares


K = 2
M = 4
ATOM_COUNT = K * M
HARD_INTERACTION = 24.0
MIN_OBJECTIVE_INTERACTION = 0.4
MAX_OBJECTIVE_INTERACTION = 3.2
HARD_FLOOR_INTERACTION = 6.5
ALLOWED_CEILING_INTERACTION = 4.0
FINAL_DETUNING = 5.0
MAX_AMPLITUDE = 8.0
PULSE_DURATION_NS = 5000
LAYOUT_SEED = 24098
EMULATED_SHOTS = 10_000
RESULT_SIGNIFICANT_DIGITS = 12


@dataclass(frozen=True)
class Candidate:
    net: int
    candidate: int
    cost: float
    nodes: frozenset[int]


@dataclass(frozen=True)
class Repair8:
    source: str
    candidates: tuple[Candidate, ...]
    conflict_pairs: frozenset[tuple[int, int]] | None = None

    def __post_init__(self):
        if len(self.candidates) != ATOM_COUNT:
            raise ValueError("k=2, m=4 repair must contain exactly eight candidates")
        pool_sizes = [sum(item.net == net for item in self.candidates)
                      for net in self.nets]
        if len(self.nets) != K or pool_sizes != [M, M]:
            raise ValueError("candidates must form two consecutive pools of four")

    @property
    def nets(self) -> tuple[int, ...]:
        return tuple(dict.fromkeys(item.net for item in self.candidates))

    def is_conflict(self, left: int, right: int) -> bool:
        pair = tuple(sorted((left, right)))
        if self.conflict_pairs is not None:
            return pair in self.conflict_pairs
        first, second = self.candidates[left], self.candidates[right]
        return first.net == second.net or bool(first.nodes & second.nodes)

    def feasible_choices(self) -> dict[str, float]:
        choices = {}
        for left, right in itertools.product(range(M), range(M, ATOM_COUNT)):
            if self.is_conflict(left, right):
                continue
            bits = bitstring((left, right))
            choices[bits] = (self.candidates[left].cost
                             + self.candidates[right].cost)
        return choices

    def optimal_choices(self) -> set[str]:
        feasible = self.feasible_choices()
        optimum = min(feasible.values())
        return {bits for bits, cost in feasible.items() if cost == optimum}


def bitstring(selected: tuple[int, ...]) -> str:
    return "".join("1" if index in selected else "0"
                   for index in range(ATOM_COUNT))


def load_repair(path: Path) -> Repair8:
    with path.open() as source:
        data = json.load(source)
    by_net = {net: [] for net in data["nets"]}
    for item in data["variables"]:
        by_net[item["net"]].append(Candidate(
            net=item["net"],
            candidate=item["candidate"],
            cost=float(item["cost"]),
            nodes=frozenset(item["nodes"]),
        ))
    selected = []
    for net in data["nets"][:K]:
        candidates = sorted(by_net[net], key=lambda item: item.candidate)
        if len(candidates) < M:
            raise ValueError(f"net {net} has fewer than four candidates")
        selected.extend(candidates[:M])
    repair = Repair8(str(path), tuple(selected))
    if len(repair.feasible_choices()) < 2:
        raise ValueError("eight-candidate repair has fewer than two feasible choices")
    if len(set(repair.feasible_choices().values())) < 2:
        raise ValueError("eight-candidate repair has no cost discrimination")
    return repair


def load_manifest_case(case: dict) -> Repair8:
    candidates = tuple(Candidate(
        net=item["net"],
        candidate=item["candidate"],
        cost=float(item["cost"]),
        nodes=frozenset(),
    ) for item in case["candidates"])
    return Repair8(
        source=case["source_file"],
        candidates=candidates,
        conflict_pairs=frozenset(tuple(pair) for pair in case["conflict_pairs"]),
    )


def target_interactions(repair: Repair8) -> dict[tuple[int, int], float]:
    feasible = repair.feasible_choices()
    low, high = min(feasible.values()), max(feasible.values())
    interactions = {}
    for left in range(ATOM_COUNT):
        for right in range(left + 1, ATOM_COUNT):
            if repair.is_conflict(left, right):
                interactions[left, right] = HARD_INTERACTION
            else:
                state = bitstring((left, right))
                fraction = ((feasible[state] - low) / (high - low)
                            if high > low else 0.0)
                interactions[left, right] = (
                    MIN_OBJECTIVE_INTERACTION
                    + fraction * (MAX_OBJECTIVE_INTERACTION
                                  - MIN_OBJECTIVE_INTERACTION)
                )
    return interactions


def fit_layout(repair: Repair8, attempts: int = 64) -> tuple[np.ndarray, float]:
    interactions = target_interactions(repair)
    pairs = tuple(interactions)
    objective_targets = np.array([
        (AnalogDevice.interaction_coeff / interactions[pair]) ** (1.0 / 6.0)
        for pair in pairs
    ])
    hard = np.array([repair.is_conflict(*pair) for pair in pairs])
    hard_max_distance = (
        AnalogDevice.interaction_coeff / HARD_FLOOR_INTERACTION
    ) ** (1.0 / 6.0)
    allowed_min_distance = (
        AnalogDevice.interaction_coeff / ALLOWED_CEILING_INTERACTION
    ) ** (1.0 / 6.0)

    def distances(flat_coordinates):
        coordinates = flat_coordinates.reshape(ATOM_COUNT, 2)
        return np.array([
            np.linalg.norm(coordinates[left] - coordinates[right])
            for left, right in pairs
        ])

    def target_residual(flat_coordinates):
        actual = distances(flat_coordinates)
        return (actual[~hard] - objective_targets[~hard]) / objective_targets[~hard]

    def constrained_residual(flat_coordinates):
        coordinates = flat_coordinates.reshape(ATOM_COUNT, 2)
        actual = distances(flat_coordinates)
        minimum_shortfall = np.maximum(
            0.0, AnalogDevice.min_atom_distance - actual
        )
        hard_excess = np.maximum(0.0, actual[hard] - hard_max_distance)
        allowed_shortfall = np.maximum(
            0.0, allowed_min_distance - actual[~hard]
        )
        radius_excess = np.maximum(
            0.0, np.linalg.norm(coordinates, axis=1)
            - AnalogDevice.max_radial_distance
        )
        return np.concatenate((
            0.25 * target_residual(flat_coordinates),
            20.0 * minimum_shortfall / AnalogDevice.min_atom_distance,
            10.0 * hard_excess / hard_max_distance,
            10.0 * allowed_shortfall / allowed_min_distance,
            10.0 * radius_excess / AnalogDevice.max_radial_distance,
        ))

    rng = np.random.default_rng(LAYOUT_SEED)
    best = None
    square = np.array([
        [-2.5, -2.5], [-2.5, 2.5], [2.5, -2.5], [2.5, 2.5]
    ])
    for attempt in range(attempts):
        initial = np.vstack((
            square + [-8.0, 0.0],
            -square + [8.0, 0.0],
        ))
        if attempt:
            initial += rng.normal(0.0, 1.0, size=initial.shape)
        result = least_squares(
            constrained_residual, initial.ravel(), bounds=(-35.0, 35.0),
            max_nfev=10_000
        )
        coordinates = result.x.reshape(ATOM_COUNT, 2)
        coordinates -= coordinates.mean(axis=0)
        pair_distances = [
            np.linalg.norm(coordinates[left] - coordinates[right])
            for left, right in pairs
        ]
        minimum_distance = min(pair_distances)
        if minimum_distance < AnalogDevice.min_atom_distance:
            coordinates *= AnalogDevice.min_atom_distance / minimum_distance
        score = float(np.sqrt(np.mean(np.square(
            constrained_residual(coordinates.ravel())
        ))))
        if best is None or score < best[1]:
            best = coordinates, score
    return best


def physical_energies(coordinates: np.ndarray) -> dict[str, float]:
    energies = {}
    for state in range(1 << ATOM_COUNT):
        bits = f"{state:0{ATOM_COUNT}b}"
        occupied = [index for index, value in enumerate(bits) if value == "1"]
        energy = -FINAL_DETUNING * len(occupied)
        for offset, left in enumerate(occupied):
            for right in occupied[offset + 1:]:
                distance = np.linalg.norm(coordinates[left] - coordinates[right])
                energy += AnalogDevice.interaction_coeff / distance ** 6
        energies[bits] = float(energy)
    return energies


def layout_metrics(repair: Repair8, coordinates: np.ndarray) -> dict[str, float]:
    distances = {}
    interactions = {}
    for left in range(ATOM_COUNT):
        for right in range(left + 1, ATOM_COUNT):
            pair = left, right
            distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
            distances[pair] = distance
            interactions[pair] = AnalogDevice.interaction_coeff / distance ** 6
    hard = [value for pair, value in interactions.items()
            if repair.is_conflict(*pair)]
    allowed = [value for pair, value in interactions.items()
               if not repair.is_conflict(*pair)]
    return {
        "minimum_atom_distance_um": min(distances.values()),
        "maximum_radial_distance_um": float(np.max(
            np.linalg.norm(coordinates, axis=1)
        )),
        "minimum_hard_interaction_rad_per_us": min(hard),
        "maximum_allowed_interaction_rad_per_us": max(allowed),
    }


def verify_ground_state(repair: Repair8, coordinates: np.ndarray) -> set[str]:
    energies = physical_energies(coordinates)
    ground_energy = min(energies.values())
    ground = {bits for bits, energy in energies.items()
              if math.isclose(energy, ground_energy, abs_tol=1e-7)}
    expected = repair.optimal_choices()
    if not ground or not ground <= expected:
        raise ValueError(
            f"physical ground state {ground} is not a subset of repair optima {expected}"
        )
    return ground


def build_sequence(coordinates: np.ndarray) -> Sequence:
    register = Register.from_coordinates(coordinates, center=True, prefix="q")
    sequence = Sequence(register, AnalogDevice)
    sequence.declare_channel("ising", "rydberg_global")
    amplitude = InterpolatedWaveform(
        PULSE_DURATION_NS, [1e-6, MAX_AMPLITUDE, MAX_AMPLITUDE, 1e-6]
    )
    detuning = InterpolatedWaveform(
        PULSE_DURATION_NS,
        [-FINAL_DETUNING, -FINAL_DETUNING, FINAL_DETUNING, FINAL_DETUNING],
    )
    sequence.add(Pulse(amplitude, detuning, 0.0), "ising")
    sequence.to_abstract_repr()
    return sequence


def decode_probabilities(probabilities: dict[str, float]) -> dict[str, float]:
    translation = str.maketrans({"g": "0", "r": "1"})
    return {state.translate(translation): probability
            for state, probability in probabilities.items()}


def emulate(sequence: Sequence) -> dict[str, float]:
    config = EmulationConfig(
        observables=(StateResult(),),
        default_evaluation_times=(1.0,),
        sampling_rate=0.05,
    )
    result = QutipBackendV2(sequence, config=config).run()
    return decode_probabilities(result.final_state.probabilities())


def sample_probabilities(probabilities: dict[str, float], source: str,
                         shots: int = EMULATED_SHOTS) -> dict[str, int]:
    states = sorted(probabilities)
    source_seed = int.from_bytes(
        hashlib.sha256(source.encode()).digest()[:8], "big"
    )
    rng = np.random.default_rng(LAYOUT_SEED ^ source_seed)
    sampled = rng.multinomial(shots, [probabilities[state] for state in states])
    return {state: int(count) for state, count in zip(states, sampled) if count}


def wilson_interval(successes: int, trials: int,
                    z_score: float = 1.959963984540054) -> tuple[float, float]:
    proportion = successes / trials
    denominator = 1.0 + z_score ** 2 / trials
    center = (proportion + z_score ** 2 / (2.0 * trials)) / denominator
    margin = (z_score / denominator
              * math.sqrt(proportion * (1.0 - proportion) / trials
                          + z_score ** 2 / (4.0 * trials ** 2)))
    return center - margin, center + margin


def canonicalize_numbers(value):
    if isinstance(value, float):
        return float(f"{value:.{RESULT_SIGNIFICANT_DIGITS}g}")
    if isinstance(value, dict):
        return {key: canonicalize_numbers(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [canonicalize_numbers(item) for item in value]
    return value


def analyze_repair(repair: Repair8) -> dict:
    coordinates, layout_error = fit_layout(repair)
    geometry = layout_metrics(repair, coordinates)
    physical_ground = verify_ground_state(repair, coordinates)
    exact_optima = repair.optimal_choices()
    sequence = build_sequence(coordinates)
    probabilities = emulate(sequence)
    feasible = repair.feasible_choices()
    sampled = sample_probabilities(probabilities, repair.source)
    feasible_probability = sum(probabilities.get(bits, 0.0) for bits in feasible)
    optimal_probability = sum(probabilities.get(bits, 0.0) for bits in exact_optima)
    feasible_samples = sum(sampled.get(bits, 0) for bits in feasible)
    optimal_samples = sum(sampled.get(bits, 0) for bits in exact_optima)
    conditional = optimal_probability / feasible_probability
    uniform = len(exact_optima) / len(feasible)
    return {
        "source": repair.source,
        "k": K,
        "m": M,
        "atom_count": ATOM_COUNT,
        "hilbert_dimension": 1 << ATOM_COUNT,
        "one_hot_selection_count": M ** K,
        "device": AnalogDevice.name,
        "pulser_profile_is_ruby_twin": False,
        "candidate_ids": [item.candidate for item in repair.candidates],
        "candidate_costs": [item.cost for item in repair.candidates],
        "coordinates_um": coordinates.round(6).tolist(),
        "layout_relative_rms_error": layout_error,
        **geometry,
        "sequence_duration_ns": sequence.get_duration(),
        "abstract_sequence_bytes": len(sequence.to_abstract_repr()),
        "feasible_choice_count": len(feasible),
        "exact_optimal_choices": sorted(exact_optima),
        "physical_ground_states": sorted(physical_ground),
        "physical_ground_degeneracy_recall": (
            len(physical_ground) / len(exact_optima)
        ),
        "feasible_probability": feasible_probability,
        "optimal_probability": optimal_probability,
        "conditional_optimal_gain": conditional / uniform,
        "emulated_shots": EMULATED_SHOTS,
        "sampling_seed": LAYOUT_SEED,
        "sampled_feasible_count": feasible_samples,
        "sampled_feasible_wilson_95": wilson_interval(
            feasible_samples, EMULATED_SHOTS
        ),
        "sampled_optimal_count": optimal_samples,
        "sampled_optimal_wilson_95": wilson_interval(
            optimal_samples, EMULATED_SHOTS
        ),
        "probabilities": dict(sorted(probabilities.items())),
    }


def analyze_manifest(path: Path) -> dict:
    with path.open() as source:
        manifest = json.load(source)
    rows = [analyze_repair(load_manifest_case(case))
            for case in manifest["cases"]]
    result = {
        "schema": "qoblib-pulser-eight-atom-results-v1",
        "manifest": str(path),
        "selection_rule": manifest["selection_rule"],
        "runtime": {
            "python": platform.python_version(),
            "pulser": importlib.metadata.version("pulser"),
            "pulser-simulation": importlib.metadata.version("pulser-simulation"),
            "numpy": importlib.metadata.version("numpy"),
            "scipy": importlib.metadata.version("scipy"),
        },
        "case_count": len(rows),
        "validated_case_count": len(rows),
        "cases": rows,
    }
    if rows:
        for field in ("feasible_probability", "optimal_probability",
                      "conditional_optimal_gain"):
            values = [row[field] for row in rows]
            result[field] = {
                "min": min(values),
                "median": statistics.median(values),
                "max": max(values),
            }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = canonicalize_numbers(analyze_manifest(args.manifest))
    text = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())