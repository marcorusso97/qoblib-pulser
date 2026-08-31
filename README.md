# QOBLIB Steiner repair on Pulser

This repository contains a reproducible eight-atom feasibility benchmark for a
neutral-atom repair subproblem derived from the QOBLIB Steiner tree packing
instances.

The benchmark keeps four candidate trees for each of two nets (`k = 2`,
`m = 4`). It fits an eight-atom register, validates the register and global
Rydberg sequence against Pulser's public `AnalogDevice` profile, enumerates the
complete final Hamiltonian spectrum, and simulates the 256-dimensional state
space with `pulser-simulation`.

Five frozen repairs are included. All five pass the geometric, spectral, and
Pulser validation checks. The median exact feasible probability is 50.94%, the
median exact optimal probability is 32.88%, and the median conditional gain over
uniform feasible sampling is 1.147x.

The full formulation and per-instance results are documented in
[FEASIBILITY.md](FEASIBILITY.md).

## Reproduce

Python 3.12 and [uv](https://docs.astral.sh/uv/) are required.

```bash
uv sync --group test
uv run pytest -q
uv run python benchmark.py cases.json --out results.json
```

The final command deterministically regenerates `results.json`. Runtime versions,
atom coordinates, validation metrics, exact probabilities, and finite-shot
confidence intervals are recorded in that file.

## Files

- `benchmark.py`: reduction, layout fitting, validation, emulation, and metrics
- `cases.json`: frozen inputs and source hashes
- `results.json`: complete emulator output
- `tests/test_benchmark.py`: focused correctness and device-validation tests
- `FEASIBILITY.md`: technical feasibility report

## Scope

The benchmark targets Pulser's public `AnalogDevice` profile. That profile has
not been confirmed as a digital twin of the RUBY system. The repository therefore
supports an SDK and ideal-emulator feasibility claim, not a hardware-performance
or quantum-advantage claim.

## Data provenance

The frozen repair data were derived from the
[Quantum Optimization Benchmarking Library](https://github.com/ZIB-AOPT/QOBLIB)
Steiner tree packing benchmark. Source filenames and SHA-256 hashes are retained
in `cases.json`.
