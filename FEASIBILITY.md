# Eight-Atom Pulser Feasibility Proof for QOBLIB Steiner Repairs

Implementation and frozen artifacts:
[github.com/marcorusso97/qoblib-pulser-feasibility](https://github.com/marcorusso97/qoblib-pulser-feasibility)

## Result

**Passed on five frozen, QOBLIB-derived `k = 2`, `m = 4` repairs under Pulser's
public `AnalogDevice` profile.**

This proof uses eight atoms, 16 one-hot candidate selections, and the complete
256-dimensional Hilbert space. All five registers and pulse schedules validate, every
physical ground state is an exact minimum-cost repair, and the same frozen global pulse
schedule produces both a measurable feasibility signal and an optimization signal.

`AnalogDevice` is a public Pasqal/Pulser profile, not a confirmed RUBY digital twin. The
result is therefore an SDK and ideal-emulator proof of concept, not yet a RUBY-matched
hardware feasibility certificate.

## Why eight atoms are non-trivial

Two original LNS nets retain all four candidates:

$$
N = k m = 2 \times 4 = 8, \qquad
|\mathcal{F}_{\mathrm{one-hot}}| = 4^{2} = 16, \qquad
\dim \mathcal{H} = 2^{8} = 256.
$$

The emulator evolves the full 256-state system. Feasibility is not imposed during
simulation or decoding: probability outside the valid one-hot, conflict-free states is
reported as infeasible.

## Frozen benchmark selection

The hardware-native cohort was defined before observing any emulator outcome:

> Select the first five repair files in lexical order whose first two four-candidate
> pools have distinct feasible costs and at most three cross-pool conflicts. Accept a
> fitted geometry only when every physical ground state is an exact repair optimum.

This produced `sub_00013`, `sub_00018`, `sub_00031`, `sub_00032`, and `sub_00035`.
The cohort contains 14-16 feasible choices per repair; `sub_00031` includes a real
cross-pool conflict. No case was selected or removed using emulator probabilities.

This deterministic screening is part of the neutral-atom co-design workflow and can be
reproduced from the tracked manifest and source hashes.

## Hamiltonian and encoding

The emulator uses

$$
H(t) = \frac{\Omega(t)}{2} \sum_i \sigma_i^x
- \delta(t) \sum_i n_i
+ \sum_{i < j} \frac{C_6}{r_{ij}^{6}} n_i n_j.
$$

Same-net pairs and cross-net conflicts are assigned to the high-interaction class.
Allowed cross-net pairs are assigned to a lower-interaction class ordered by their
classical candidate cost. The deterministic layout fitter first enforces the energetic
separation

$$
V_{\mathrm{hard}} > \delta_f > V_{\mathrm{allowed}},
$$

and only then minimizes cost-order distortion. The shared schedule uses a 5,000 ns
global Rydberg pulse, maximum amplitude `8 rad/us`, and detuning from `-5` to
`+5 rad/us`. Layout seed and sampling seed are `24098`.

For the five accepted layouts:

- minimum atom distance: `5.000 um`;
- maximum radial distance: `9.421 um`, below the 38 um profile limit;
- minimum hard interaction: `6.479-6.481 rad/us`;
- maximum allowed interaction: `4.006-4.016 rad/us`;
- final detuning: `5 rad/us`;
- schedule duration: `5,000 ns`, below the 6,000 ns profile limit.

All 256 final classical energies are enumerated before dynamic emulation. Geometry-induced
breaking of exact degeneracy is allowed only when the surviving physical ground states
are a non-empty subset of the exact QOBLIB optima. Any infeasible or suboptimal physical
ground state rejects the case.

## Emulator results

Exact probabilities come from the final statevector. Finite-shot values use 10,000
deterministic emulator samples per case; intervals are two-sided Wilson 95% intervals.
The gain is the probability of an exact optimum conditioned on feasibility, divided by
uniform sampling over the feasible choices.

| Repair | Feasible choices | Exact feasible | Sampled feasible (95% CI) | Exact optimum | Sampled optimum (95% CI) | Conditional gain vs uniform |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `sub_00013` | 16 | 49.55% | 49.82% (48.84-50.80%) | 21.13% | 21.30% (20.51-22.11%) | 1.137x |
| `sub_00018` | 16 | 50.94% | 50.50% (49.52-51.48%) | 32.88% | 33.00% (32.09-33.93%) | 1.148x |
| `sub_00031` | 15 | 57.81% | 57.95% (56.98-58.91%) | 41.48% | 41.05% (40.09-42.02%) | 1.196x |
| `sub_00032` | 16 | 50.94% | 50.98% (50.00-51.96%) | 32.88% | 32.49% (31.58-33.41%) | 1.147x |
| `sub_00035` | 16 | 50.94% | 50.68% (49.70-51.66%) | 32.88% | 32.57% (31.66-33.50%) | 1.147x |
| **Median** | **16** | **50.94%** | - | **32.88%** | - | **1.147x** |

The optimization gain exceeds 1 in all 5/5 cases. This is an ideal-emulator signal on a
small, screened corpus, not evidence of quantum advantage.

## Reproduction

```bash
uv sync --group test
uv run pytest -q
uv run python benchmark.py cases.json --out results.json
```

Tracked artifacts:

- [`benchmark.py`](benchmark.py);
- [`cases.json`](cases.json);
- [`results.json`](results.json);
- [`tests/test_benchmark.py`](tests/test_benchmark.py).

Runtime recorded in the result artifact: Python 3.12.13, Pulser 1.9.0,
`pulser-simulation` 1.9.0, NumPy 2.5.2, and SciPy 1.18.1.

## Demonstrated scope

The result establishes a non-trivial eight-atom proof of concept for geometry-screened
QOBLIB repairs under global Rydberg control. The complete workflow succeeds for all five
frozen cases: reduction, exact enumeration, geometry fitting, device validation,
serialization, full-state emulation, deterministic decoding, and finite-shot reporting.

Before calling this a RUBY feasibility proof, the hosting entity must confirm RUBY's
accepted Pulser device/control model. The next scientific gate is a preregistered
12-atom geometry-native cohort under that hardware-matched profile, including noise and
atom-loss models supplied or approved by the provider.