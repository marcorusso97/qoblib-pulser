import importlib.util
import json
from pathlib import Path

import numpy as np


PATH = Path(__file__).parents[1] / "benchmark.py"
SPEC = importlib.util.spec_from_file_location("benchmark", PATH)
PROOF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PROOF)


def write_repair(path):
    variables = []
    for net, base_cost, node_offset in ((2, 10, 0), (5, 20, 100)):
        for candidate in range(4):
            variables.append({
                "net": net,
                "candidate": candidate,
                "cost": base_cost + candidate,
                "nodes": [node_offset + candidate],
            })
    path.write_text(json.dumps({"nets": [2, 5], "variables": variables}))


def test_reduction_has_sixteen_weighted_choices(tmp_path):
    path = tmp_path / "repair.json"
    write_repair(path)

    repair = PROOF.load_repair(path)

    assert repair.nets == (2, 5)
    assert len(repair.feasible_choices()) == 16
    assert repair.optimal_choices() == {"10001000"}
    assert repair.feasible_choices()["00010001"] == 36


def test_eight_atom_layout_sequence_and_spectrum(tmp_path):
    path = tmp_path / "repair.json"
    write_repair(path)
    repair = PROOF.load_repair(path)

    coordinates, error = PROOF.fit_layout(repair, attempts=8)
    metrics = PROOF.layout_metrics(repair, coordinates)

    assert coordinates.shape == (8, 2)
    assert error < 0.25
    assert metrics["minimum_atom_distance_um"] >= 5.0
    assert metrics["maximum_radial_distance_um"] <= 38.0
    assert metrics["minimum_hard_interaction_rad_per_us"] >= 6.4
    assert metrics["maximum_allowed_interaction_rad_per_us"] <= 4.1
    assert PROOF.verify_ground_state(repair, coordinates) == {"10001000"}
    sequence = PROOF.build_sequence(coordinates)
    assert sequence.get_duration() == PROOF.PULSE_DURATION_NS
    abstract = json.loads(sequence.to_abstract_repr())
    assert len(abstract["register"]) == 8
    assert abstract["device"]["name"] == "AnalogDevice"


def test_ground_state_gate_rejects_nonoptimal_state(monkeypatch, tmp_path):
    path = tmp_path / "repair.json"
    write_repair(path)
    repair = PROOF.load_repair(path)
    energies = {f"{state:08b}": 1.0 for state in range(256)}
    energies["01000100"] = 0.0
    monkeypatch.setattr(PROOF, "physical_energies", lambda _coordinates: energies)

    try:
        PROOF.verify_ground_state(repair, None)
    except ValueError as error:
        assert "not a subset" in str(error)
    else:
        raise AssertionError("nonoptimal physical ground state was accepted")


def test_exact_optimum_probability_keeps_degenerate_optima(monkeypatch, tmp_path):
    path = tmp_path / "repair.json"
    write_repair(path)
    data = json.loads(path.read_text())
    data["variables"][1]["cost"] = 10
    path.write_text(json.dumps(data))
    repair = PROOF.load_repair(path)
    probabilities = {"10001000": 0.3, "01001000": 0.2, "00000000": 0.5}
    class FakeSequence:
        def get_duration(self):
            return 5000

        def to_abstract_repr(self):
            return "{}"

    monkeypatch.setattr(
        PROOF, "fit_layout", lambda _repair: (np.zeros((8, 2)), 0.0)
    )
    monkeypatch.setattr(PROOF, "layout_metrics", lambda *_args: {})
    monkeypatch.setattr(
        PROOF, "verify_ground_state", lambda _repair, _coordinates: {"10001000"}
    )
    monkeypatch.setattr(PROOF, "build_sequence", lambda _coordinates: FakeSequence())
    monkeypatch.setattr(PROOF, "emulate", lambda _sequence: probabilities)

    row = PROOF.analyze_repair(repair)

    assert row["physical_ground_states"] == ["10001000"]
    assert row["optimal_probability"] == 0.5


def test_sampling_is_reproducible_and_wilson_contains_estimate():
    probabilities = {"00000000": 0.4, "10001000": 0.6}

    first = PROOF.sample_probabilities(probabilities, "case", shots=1000)
    second = PROOF.sample_probabilities(probabilities, "case", shots=1000)
    interval = PROOF.wilson_interval(first["10001000"], 1000)

    assert first == second
    assert sum(first.values()) == 1000
    assert interval[0] < first["10001000"] / 1000 < interval[1]