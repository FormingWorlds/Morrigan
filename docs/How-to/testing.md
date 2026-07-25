# Testing

Morrigan's suite follows the PROTEUS ecosystem test standards: tiered markers, physics-invariant assertions on every physics routine, and reference-pinned anchors recorded per source file under [Validation anchors](../Validation/index.md).

[![tests](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2FFormingWorlds%2FMorrigan%2Fbadges%2Ftests-total.json)](https://github.com/FormingWorlds/Morrigan/tree/badges)
[![unit tests](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2FFormingWorlds%2FMorrigan%2Fbadges%2Ftests-unit.json)](https://github.com/FormingWorlds/Morrigan/tree/badges)
[![integration tests](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2FFormingWorlds%2FMorrigan%2Fbadges%2Ftests-integration.json)](https://github.com/FormingWorlds/Morrigan/tree/badges)

These counts are regenerated from the suite on every push to `main`. They report how many tests exist, not whether they pass; the pass or fail state is the workflow status badge on the [front page](../index.md). The "integration tests" count covers the `smoke`, `integration` and `slow` tiers together, so it does not match any single row of the table below: most of it is the `smoke` tier, which runs on every pull request rather than nightly. Counting is by collection, so a test that is collected but skipped at run time, through `importorskip` or `skipif`, still counts.

## Tiers

| Marker | What it tests | Budget | When CI runs it |
|---|---|---|---|
| `unit` | Single functions, small systems | < 100 ms each | Every pull request |
| `smoke` | Real full-system runs, one small system | < 30 s each | Every pull request |
| `integration` | Seed ensembles, statistics over many runs | Minutes | Nightly |
| `slow` | Full physics validation | Up to hours | Nightly |

Every test file carries a module-level `pytestmark` naming its tier and timeout; tests without a tier marker are invisible to CI.

## Commands

```bash
pytest                                          # everything
pytest -m "(unit or smoke) and not skip"        # the fast tiers (what PRs run)
pytest -m "integration and not skip"            # the nightly ensembles
pytest --cov=morrigan                           # with coverage
```

Structure and quality checks, both enforced on pull requests:

```bash
bash tools/validate_test_structure.sh           # marker validation
python tools/check_test_quality.py --check      # anti-happy-path lint
```

## Coverage gates

Both gates are one-way ratchets capped at the 90 percent ecosystem ceiling, updated with `python tools/update_coverage_threshold.py` and never manually decreased:

- the fast gate (`[tool.morrigan.coverage_fast]`) is checked against the unit and smoke tiers on every pull request;
- the full gate (`[tool.coverage.report]`) is checked against the complete suite by the nightly workflow, which also uploads coverage to Codecov.

## Writing new tests

Tests mirror source one to one: `src/morrigan/<file>.py` is exercised by `tests/test_<file>.py`. One documented exception exists: the cross-cutting ensemble statistics (`tests/test_ensemble_statistics.py`). Every test on a physics routine asserts at least one physical invariant (conservation, positivity or boundedness, monotonicity or symmetry, or a pinned value with a discrimination guard), every test carries a docstring naming the scenario it verifies, and any test that triggers a random draw seeds the global numpy state explicitly and states the seed. Optional dependencies are import-or-skipped at module top.
