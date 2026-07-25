# Morrigan AI Agent Guidelines

**Trust these instructions.** Only search if information is incomplete or found to be in error.

**Identity & Mission**: You are an expert Scientific Software Engineer working on the Morrigan module of the PROTEUS ecosystem.

## High-Level Instructions

> ### Rule files you MUST read on every session
>
> Morrigan keeps its Claude-Code rule files under `.github/.claude/rules/` (NOT the conventional repo-root `.claude/`, which is gitignored and so cannot be shared with collaborators). Claude Code does NOT auto-discover the rules at this unusual path. Read them explicitly at the start of every session and any time you open a related file:
>
> - [`.github/.claude/rules/morrigan-tests.md`](.claude/rules/morrigan-tests.md) -- test quality deep-dive: anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, global-RNG seeding discipline, timeline-chain invariants. **Required reading before editing any file under `tests/**` or `src/morrigan/**`.**
> - [`.github/.claude/rules/morrigan-code-review.md`](.claude/rules/morrigan-code-review.md) -- review-pass gate, domain-aware physics review (orbital element bounds, merger mass closure, unit boundaries, PROTEUS record-schema contract). **Required reading before any code review pass.**
>
> These two files plus this one are the canonical sources of truth for testing rigor and review criteria.

1. **Always** read the two rule files above plus the Testing Standards section below before any code change.
2. **Always** inform the user that you are reading in this file by printing a message at the start of your response: "(Read in copilot-instructions.md...)"
3. When creating a PR, **always** state what changed and why in terms of the model, with a Testing section describing how it was verified.
4. **Claude-specific**: `CLAUDE.md` is a symlink to this file. Session learnings, plans, and memories live in `~/.claude/projects/<repo>/memory/`; they do NOT live in this repository.

## Ecosystem Context

Morrigan is the giant-impact accretion module of the PROTEUS ecosystem: a semi-analytical Monte Carlo model for the post-disk dynamical evolution of planetary systems, following Kimura et al. (2025). It evolves a system of planetary embryos through orbital crossings, gravitational scatterings, and giant impacts, and reports the impact history of a selected survivor. The main [PROTEUS](https://github.com/FormingWorlds/PROTEUS) coupled atmosphere-interior framework consumes that history through `morrigan.run_system` to drive its accretion coupling: each reported impact grows the PROTEUS planet, re-melts its mantle, delivers volatiles, strips atmosphere, and moves its orbit.

Sister modules: AGNI (radiative transfer), SOCRATES (spectral radiative transfer), JANUS (1D atmosphere), MORS (stellar evolution), CALLIOPE (outgassing), ARAGOG / SPIDER (interior), VULCAN (chemistry), ZEPHYRUS / BOREAS (escape), Zalmoxis (structure), Obliqua (tides).

**Project Type**: Scientific simulation module (Python).

**Languages**: Python 3.11+.

**Size**: 11 source files in `src/morrigan/`, ~1.3k LOC.

**Target Runtime**: Python 3.11+ on Linux / macOS.

## Build & Validation

### Environment Setup

**Developer Install**:

```bash
git clone git@github.com:FormingWorlds/Morrigan.git
cd Morrigan
pip install -e ".[develop]"
pre-commit install -f
```

Morrigan has no compiled dependencies. The `pylaplace` dependency supplies the Laplace-Lagrange secular machinery and installs from PyPI.

### Test Commands

**Run all tests**:

```bash
pytest
```

**Run by category** (matches CI):

```bash
pytest -m "(unit or smoke) and not skip and not slow and not integration"   # What PR checks run
pytest -m unit                    # Fast unit tests (< 100 ms each)
pytest -m integration             # Ensemble statistics over many seeds (nightly)
pytest -m slow                    # Full physics validation (nightly)
```

**With coverage**:

```bash
pytest -m "(unit or smoke) and not skip" --cov=morrigan
coverage report
```

**Coverage thresholds** (in `pyproject.toml`; a one-way ratchet raised with `tools/update_coverage_threshold.py` when coverage improves, never manually decreased, capped at 90):

- Fast gate (`[tool.morrigan.coverage_fast]`, unit + smoke, every PR): ratcheting toward **90%** (the PROTEUS-ecosystem ceiling).
- Full gate (`[tool.coverage.report]`, unit + smoke + integration + slow, nightly): ratcheting toward **90%**.

**Validate test structure**:

```bash
bash tools/validate_test_structure.sh
```

**Test quality lint** (blocking on PRs):

```bash
python tools/check_test_quality.py --check
```

### Lint Commands

**Always run before committing**:

```bash
ruff check src/ tests/ tools/
ruff check --fix src/ tests/ tools/
ruff format src/ tests/ tools/
```

### Validation Pipeline

**CI runs on PRs** (`.github/workflows/tests.yaml`):

1. **Unit + smoke tests**: `pytest -m "(unit or smoke) and not skip and not slow and not integration" --cov=morrigan`.
2. **Fast coverage gate**: `[tool.morrigan.coverage_fast].fail_under` checked against the unit + smoke coverage.
3. **Test structure**: `bash tools/validate_test_structure.sh`.
4. **Test quality**: `python tools/check_test_quality.py --check` (blocking).
5. **Coverage ratchet guard**: rejects any PR that lowers `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.
6. **Lint**: `ruff check src/ tests/ tools/`.

The test job carries `timeout-minutes: 10`; a PR tier that cannot finish inside it belongs in the nightly tier.

**Nightly CI** (`.github/workflows/nightly.yml`):

- Slow tiers: `pytest -m "(integration or slow) and not skip"` (the seed-ensemble statistics live here).
- Triggered nightly at 03:00 UTC plus `workflow_dispatch` for manual runs.

## Project Layout

### Key Directories

- `src/morrigan/` - Main Python source code (flat layout, 12 files)
  - `__init__.py` - Public API re-exports (`run_system`) (utility)
  - `_version.py` - Auto-generated by setuptools-scm (utility)
  - `constants.py` - Physical constants and unit conversions (utility)
  - `sort_planet.py` - Dead-planet removal and semi-major-axis sorting (utility)
  - `driver.py` - Orchestration, config handling, per-impact record building, survivor selection (physics)
  - `orbit_cross_K25.py` - Orbit-crossing Monte Carlo event loop (physics)
  - `interaction_timescales.py` - Collision and viscous-stirring timescales (physics)
  - `crossing_pair.py` - Interacting-pair identification and crossing times (physics)
  - `merge_embryo.py` - Merger bookkeeping: masses, orbits, eccentricities (physics)
  - `secular_solution.py` - Laplace-Lagrange secular eigenmodes via pylaplace (physics)
  - `helper_functions.py` - Kepler period, Hill radius, escape eccentricity, Rayleigh draws (physics)

- `tests/` - Test suite. Each physics source has a 1:1 test file at `tests/test_<file>.py`. One documented exception: cross-cutting ensemble statistics (`test_ensemble_statistics.py`).

- `tools/` - Build / utility scripts
  - `check_test_quality.py` - AST linter (blocking on PRs)
  - `update_coverage_threshold.py` - One-way coverage ratchet (capped at 90)
  - `check_file_sizes.sh` - Line-cap hook on this file
  - `validate_test_structure.sh` - Module-level marker validator
  - `generate_test_badges.py` - Test-count badge JSON for the docs site and the ecosystem dashboard

- `docs/` - Documentation (Zensical; Diátaxis structure)
  - `Explanations/` - Concept pages (model, limitations, PROTEUS coupling)
  - `How-to/` - Task guides (installation, running, tests)
  - `Reference/` - Parameters and API
  - `Validation/<file>.md` - Per-source-file inventory of `@pytest.mark.reference_pinned` tests

### Configuration Files

- `pyproject.toml` - Package metadata, pytest config, coverage thresholds, ruff rules.
- `mkdocs.yml` - Documentation configuration (used by Zensical).
- `.github/workflows/` - CI / CD pipelines (`tests.yaml` PR gate, `nightly.yml` slow tiers, `docs.yaml` docs deploy, `publish-test-badges.yml` test-count badges).

### Entry Points

- **Python API**: `from morrigan import run_system`. Takes a configuration dictionary plus a seed and returns the full system history and the per-impact record list for a selected planet.
- **CLI**: `morrigan` (defined in `driver.py`), driven by a TOML settings file.

## Testing Standards

Morrigan is scientific simulation code, so the test suite is held to physics-grade rigor. The rules below are the contract; the deep-dive lives in [`.github/.claude/rules/morrigan-tests.md`](.claude/rules/morrigan-tests.md). Read that file before editing any test file or any source file under `src/morrigan/**`. The two files must be kept in sync; if you change one, mirror the change in the other.

### Structure

- Tests mirror source 1:1: `src/morrigan/<file>.py` -> `tests/test_<file>.py`. One documented exception: the cross-cutting ensemble-statistics file (`tests/test_ensemble_statistics.py`).
- Framework: `pytest` exclusively in the `tests/` directory.

### Markers and the module-level marker rule

| Marker | What it tests | Speed budget | When CI runs it |
|---|---|---|---|
| `@pytest.mark.unit` | Python logic, single functions, small systems | < 100 ms per test | Every PR |
| `@pytest.mark.smoke` | Real full-system run, one small system | < 30 s per test | Every PR |
| `@pytest.mark.integration` | Seed ensembles, statistics over many runs | Minutes per test | Nightly only |
| `@pytest.mark.slow` | Full physics validation | Up to hours per test | Nightly only |
| `@pytest.mark.skip` | Placeholder, deliberately disabled | n/a | Never |

**Mandatory module-level marker** (no exceptions): every test file begins with

```python
pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

with timeouts: 30 s for unit, 60 s for smoke, 300 s for integration, 3600 s for slow. Per-function markers are additive but do not replace the module-level marker. Tests without a tier marker are invisible to CI.

### Physics validity

Every unit test on a **physics source** (`driver.py`, `orbit_cross_K25.py`, `interaction_timescales.py`, `crossing_pair.py`, `merge_embryo.py`, `secular_solution.py`, `helper_functions.py`) must assert at least one of:

- **Conservation**: merger mass closure (`M_merged = M_target + M_impactor`, exactly, since merging is perfect), chain continuity along one body's impact history.
- **Positivity / boundedness**: masses, radii, densities, and semi-major axes strictly positive; eccentricities in `[0, 1)`; impact parameter `b = sin(beta)` in `[0, 1]`; collision speed at or above the mutual escape speed.
- **Monotonicity or symmetry**: timeline strictly ordered in time; Kepler period increasing with semi-major axis; Hill radius increasing with planet mass; timescales scaling as the closed forms prescribe.
- **Pinned numeric value with a discrimination guard**: a closed-form or published value pinned via `pytest.approx`, with explicit assertions that wrong-formula results (wrong exponent, wrong mass denominator, wrong base of the exponential) land outside the tolerance.

Utility sources (`__init__.py`, `_version.py`, `constants.py`, `sort_planet.py`) are **exempt** from the physics-invariant requirement but still subject to the anti-happy-path rules.

Tag every test that asserts a physical invariant with `@pytest.mark.physics_invariant`. Per-source-file granularity: each of the seven physics files needs at least one such test in `tests/test_<file>.py`.

### Reference-pinned validation

Tag tests that pin against a published benchmark, an analytical limit, or a cross-implementation cross-check with `@pytest.mark.reference_pinned`. Each of the seven physics files must have at least one such test, recorded in `docs/Validation/<file>.md`. The `--reference-pinned-status` mode of the linter reports the punch list.

### Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case** (boundary value, two-planet minimal system, extreme mass ratio).
2. **At least one path that exercises the error contract** (documented exception, guard return, graceful clamp). If the function has no validation, exercise the limit-input behavior and assert the mathematical invariant.
3. **Assertion values that are NOT trivially derivable from the implementation**: discriminating numeric pins or property-based assertions.

**Forbidden patterns** (flagged by `tools/check_test_quality.py`): single-assert tests, standalone weak assertions, missing function docstrings, `==` adjacent to float literals, tests asserting a fixture's implicit default.

### Determinism and the global random state

Morrigan's Monte Carlo draws use numpy's **global** random state, seeded once per run by `run_system(seed=...)`. Test discipline:

- Every test that triggers a random draw seeds explicitly (`np.random.seed(42)` or a seeded `run_system` call) and states the seed in the docstring.
- Two `run_system` calls with the same seed and config must produce identical timelines; a determinism test pins this.
- Never assert on a specific random draw's value without the seed stated; never rely on test execution order for random-state setup.

### Optional-dependency imports

Any test that imports an optional dependency MUST call `pytest.importorskip('<dep>')` at module top. A dependency-light CI image will otherwise fail to collect.

### Float and numerical comparison

NEVER use `==` for floats. Use `pytest.approx(val, rel=1e-5)` or `np.testing.assert_allclose(...)`. For pinned numeric values, include a **discrimination guard**. See `morrigan-tests.md` Section 2 for the canonical pattern.

### Documentation per test

- File-level docstring: name the source under test and the invariants exercised.
- Function-level docstring: state the physical scenario or contract clause being verified. Required (lint-enforced).
- Inline comments: explain why a specific input range was chosen.

### Independent review trigger

A pull request that adds or substantially modifies > 50 lines of test code across all its commits triggers an independent review pass before merge. The denominator is PR-level (`git diff origin/main...HEAD -- 'tests/**'`).

### Tooling

- Validate test structure: `bash tools/validate_test_structure.sh`
- Test-quality lint: `python tools/check_test_quality.py --check`
- Baseline regeneration (after a deliberate sweep): `python tools/check_test_quality.py --baseline`
- Reference-pinned punch list: `python tools/check_test_quality.py --reference-pinned-status`
- Coverage ratchet (one-way, capped at 90): `python tools/update_coverage_threshold.py`

### Coverage architecture

| Gate | Tests included | Target | Enforced |
|---|---|---|---|
| Fast gate (`tool.morrigan.coverage_fast.fail_under`) | unit + smoke | Ratcheting toward **90%** | Every PR |
| Full gate (`tool.coverage.report.fail_under`) | unit + smoke + integration + slow | Ratcheting toward **90%** | Nightly |

Both gates ratchet toward 90, capped at 90 (`tools/update_coverage_threshold.py` enforces `ECOSYSTEM_CEILING = 90.0`); neither may be manually decreased.

## PROTEUS coupling contract

PROTEUS consumes `run_system` through `proteus.accretion.morrigan`. The per-impact record schema is the interface:

`time [yr], M_target_before [kg], M_impactor [kg], M_merged_after [kg], v_impact [m/s], v_esc [m/s], impact_parameter, R_target_before [m], R_impactor [m], rho_target [kg/m3], rho_impactor [kg/m3], a_before [m], a_after [m], e_after, id_target, id_impactor`

Conventions the schema guarantees:

- `v_impact` is the speed at first contact (`sqrt(v_inf^2 + v_esc^2)`), the convention the Kegerreis et al. (2020) erosion law expects.
- `M_merged_after` is the perfect-merger sum; atmospheric loss between impacts is left to the caller and bounded by PROTEUS's own timeline validation.
- Bodies carry no atmosphere at all, in the record or in the model; masses and radii are bulk values.
- Densities are recovered from mass and radius, so mass, radius, and density in one record are exactly self-consistent.
- Successive records for one body chain: each `M_target_before` equals the previous `M_merged_after` exactly, since nothing removes mass between impacts.

Any change to the record fields, units, or conventions is a breaking interface change for PROTEUS and must be flagged in the PR description.

## Code Quality

**Style** (enforced by ruff): line length < 96, `snake_case` functions, `UPPER_CASE` constants, standard type hints, docstrings with brief physical descriptions.

**Pre-commit**: runs `ruff check --fix` and the file-size hook automatically.

## Common Workflows

### Making a Code Change

1. **Create branch**: `git checkout -b <initials>/<short-description>`.
2. **Make changes** in `src/morrigan/`.
3. **Write / update tests** in `tests/test_<file>.py` (mirror structure).
4. **Run tests locally**: `pytest -m "(unit or smoke) and not skip"`.
5. **Lint**: `ruff check --fix src/ tests/ tools/`.
6. **Validate structure**: `bash tools/validate_test_structure.sh`.
7. **Test quality**: `python tools/check_test_quality.py --check`.
8. **Commit**: plain-language subject, first-person voice.
9. **Push**: CI runs automatically on PR.

### Debugging Test Failures

```bash
pytest -v --showlocals
pytest -x
pytest tests/test_<file>.py::test_function
pytest --pdb
```

## Documentation References

- **Testing rules**: `.github/.claude/rules/morrigan-tests.md`, `.github/.claude/rules/morrigan-code-review.md`
- **Docs site**: `docs/` (build and serve with `zensical serve`)
- **Model**: Kimura et al. (2025), ApJ 989, 109 (Paper I) and the application paper (Paper II)

## Project memory and session learnings

Session-specific knowledge lives outside this repository, in the Claude memory tree under `~/.claude/projects/<repo>/memory/`. What lives in this repository: this file, the two rule files, PR descriptions, commit messages, and `docs/Validation/<file>.md` pages. Do not introduce a new in-repo memory or decisions-log file.

---

## Quick Reference

```bash
# Setup
pip install -e ".[develop]"
pre-commit install -f

# Test
pytest -m "(unit or smoke) and not skip"
pytest --cov=morrigan --cov-report=html

# Lint
ruff check --fix src/ tests/ tools/
ruff format src/ tests/ tools/

# Validate
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check

# Serve docs locally
pip install -e '.[docs]'
zensical serve
```

**Remember**: Trust these instructions. Only search if information is incomplete or found to be in error.

---

> **⚠️ FILE SIZE LIMIT: This file must stay below 500 lines.** Enforced by pre-commit hook (`tools/check_file_sizes.sh`). File located at `.github/copilot-instructions.md`.
