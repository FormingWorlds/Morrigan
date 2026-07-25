---
description: Morrigan test quality deep-dive. Anti-happy-path patterns, discriminating-value guards, physics-invariant tiering, validation certification markers, global-RNG seeding discipline, timeline-chain invariants. Extends the Testing Standards section in `.github/copilot-instructions.md`.
---

# Morrigan Test Quality Rules

This file is the canonical deep-dive on test quality. The high-level summary lives in [`.github/copilot-instructions.md`](../../copilot-instructions.md) under "Testing Standards". The two files MUST stay in sync. If you change one, mirror the change in the other.

> **Discovery note.** Morrigan keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `morrigan-code-review.md` explicitly. **When opening or editing any file under `tests/**` or `src/morrigan/**`, read this file first.**

Morrigan is scientific simulation code and the test suite is held to physics-grade rigor. A test that asserts the wrong thing, or that passes for the wrong reason, is worse than no test because it generates false confidence.

---

## 1. Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case**: a boundary value (`b = 0` head-on, `b = 1` grazing, `e = 0` circular, a two-planet minimal system, an extreme mass ratio), an empty input, or an extreme physical parameter.
2. **At least one path that exercises the error contract**:
   - If the function under test has documented validation, test that the error fires AND that no side effect ran.
   - If the function is closed-form mathematics (Kepler period, Hill radius, timescales), exercise the **limit-input behavior** (equal masses, coincident orbits, zero eccentricity) and assert the corresponding mathematical invariant.
   - "No validation in source therefore no error test" is not an exemption; the limit-input substitute is.
3. **Assertion values NOT trivially derivable from the implementation**: discriminating numeric pins (Section 2) or property-based assertions (monotonicity, conservation, symmetry, boundedness).

### Forbidden patterns

These are flagged by `tools/check_test_quality.py` and rejected at PR time.

- **Single-assert test functions.** Exception: a single assertion of a hard-fail invariant (mass closure within `1e-12`) is acceptable if the test is the only test of that invariant in the file.
- **Weak assertions standing alone as the sole meaningful check**: `assert result is not None`, `assert result > 0`, `assert len(result) > 0`, `assert isinstance(result, dict)`, `assert result is None` on an implicit-None return. The three-class discrimination guard's secondary lines (sign guard, scale guard) are NOT flagged when paired with a stronger primary assertion in the same test.
- **Tests with no function-level docstring.**
- **`==` adjacent to a float literal.** Use `pytest.approx` or `np.testing.assert_allclose`.
- **Tests asserting on a fixture's implicit default.**

---

## 2. Discriminating test values

The test contract: a regression that introduces a plausible bug must fail the test. Plausible bugs in this codebase are known from its history: a wrong mass denominator in a ratio, `np.exp` where `10**` is required, a wrong pair index in a mass product, an unclamped fitted power law, a mass assigned into a fraction array. Pick input values where the wrong-formula result is far from the correct one.

### Bad / good examples

| Pattern | Bad (any-formula-passes) | Good (discriminates) |
|---|---|---|
| Crossing timescale with `10**` fit | Inputs where `exp` and `10**` nearly coincide | Inputs where the two bases differ by orders of magnitude |
| Kepler period | `a` chosen so numeric coincidences mask a wrong exponent | Two semi-major axes spanning a decade; assert the `a^(3/2)` scaling between them |
| Merger bookkeeping | Equal masses (mass-index errors degenerate) | Unequal masses on the two bodies |

### Discrimination guard (REQUIRED for pinned-value tests)

When a test pins a numeric value, include explicit assertions that the wrong-formula result would differ for each plausible bug class:

1. **Exponent or factor error**: `abs(val - wrong_value)` outside tolerance.
2. **Sign error**: assert the sign explicitly.
3. **Unit-conversion error**: pin the absolute scale with the unit named in a comment (seconds vs years and metres vs AU are the recurring boundaries here).
4. **Wrong-denominator / wrong-index selection**: when a ratio or a pair indexing is the formula's core, the guard MUST include the value the plausible wrong selection would produce.

**Argument-order pins.** A difference-based guard (`abs(f(a) - f(b)) > delta`) does NOT pin argument order: a dispatch that permutes its arguments permutes the two values and leaves the difference unchanged. Pin each side to its absolute value instead, plus the direction.

**Carve-out for conservation-style invariants.** When the primary assertion IS a conservation closure, the exponent guard is satisfied by the closure itself; sign and scale guards remain mandatory.

---

## 3. Physics-invariant assertions (tiered)

### When required

Every unit test on a **physics source** must assert at least one of the four invariant families. Physics sources:

```
src/morrigan/driver.py
src/morrigan/orbit_cross_K25.py
src/morrigan/interaction_timescales.py
src/morrigan/crossing_pair.py
src/morrigan/merge_embryo.py
src/morrigan/secular_solution.py
src/morrigan/helper_functions.py
```

Utility sources are exempt from the physics-invariant requirement but still subject to all anti-happy-path rules:

```
src/morrigan/__init__.py       (re-exports)
src/morrigan/_version.py       (auto-generated by setuptools-scm)
src/morrigan/constants.py      (pure constants, no derivation)
src/morrigan/sort_planet.py    (dead-planet removal and sorting, pure bookkeeping)
```

### The four invariant families, in Morrigan terms

1. **Conservation**
   - Merger mass closure: `M_merged_after = M_target_before + M_impactor`, exactly. Merging is perfect and the model sheds nothing, so this is an equality, not a bound.
   - Chain continuity: successive impacts on one body satisfy `M_target_before(n+1) == M_merged_after(n)` exactly; a body only ever grows, and nothing removes mass between impacts.
2. **Positivity / boundedness**
   - Masses, radii, densities, semi-major axes strictly positive; eccentricities in `[0, 1)`.
   - `b = sin(beta)` in `[0, 1]`.
   - `v_impact >= v_esc` (the contact speed carries the mutual escape speed as a floor).
   - Timeline times strictly increasing; system planet count non-increasing.
3. **Monotonicity or symmetry**
   - Kepler period scaling `P propto a^(3/2)` at fixed stellar mass.
   - Hill radius increasing with planet mass and semi-major axis.
   - Escape eccentricity symmetric under swapping the two bodies.
4. **Pinned numeric value with a discrimination guard**: Section 2. Acceptable as the sole invariant when a closed form or a published value is the contract.

### Validation certification markers

- **`@pytest.mark.physics_invariant`**: the test asserts at least one of the four families. Tag every qualifying test in a physics-source test file.
- **`@pytest.mark.reference_pinned`**: the test pins behavior against a **published benchmark** (Kimura et al. 2025 statistics and closed forms), an **analytical limit** (Kepler's third law, the mutual-Hill-radius closed form, Laplace-Lagrange two-planet eigenfrequencies), or a **cross-implementation cross-check**.
  - **Per-source-file**: each of the seven physics files needs at least one `reference_pinned` test in `tests/test_<file>.py`, recorded in `docs/Validation/<file>.md`.
  - **Status report**: `python tools/check_test_quality.py --reference-pinned-status` prints the punch list.

---

## 4. Determinism and the global random state

Morrigan's Monte Carlo layer (`rayleigh` draws, crossing-time sampling, scattering outcomes) uses numpy's **global** random state, seeded once by `run_system(seed=...)`. This makes seeding discipline the single most important determinism rule in the repo:

- **Every test that triggers any random draw seeds explicitly** and states the seed in its docstring. Either `np.random.seed(<n>)` immediately before the call, or a seeded `run_system(...)`.
- **A determinism pin exists and must stay**: two `run_system` calls with identical config and seed produce identical timelines, field for field. Any new stochastic code path must keep this test green.
- **Never assert on an unseeded draw.** Never rely on test ordering for random state: pytest may reorder, parallelise, or subset.
- **Ensemble statistics** (integration tier) sweep an explicit, fixed list of seeds, so the ensemble is reproducible and its size is visible in the test body.
- Library code must not reseed the global state outside `run_system`'s entry (a hidden `np.random.seed` in a helper would silently correlate ensemble members; this is also a review-gate item in `morrigan-code-review.md`).

---

## 5. Mocking discipline

- Default to `unittest.mock` for external calls in unit tests. Morrigan has few: file I/O on the TOML settings path is the main one.
- The Monte Carlo internals are NOT mocked in unit tests; they are seeded. A seeded draw is deterministic and real, which beats a mock that hides distributional bugs.
- A mocked function must return physically plausible values.
- NEVER mock the function under test.
- Smoke and integration tiers run the real full system.

---

## 6. Optional-dependency imports

Any test importing an optional dependency MUST call `pytest.importorskip` at module top:

```python
import pytest

hypothesis = pytest.importorskip('hypothesis')
pytest.importorskip('zephyrus.collision')
```

No test currently imports an optional dependency, and neither `hypothesis` nor
`zephyrus` is installed by any extra. The linter still recognises both in its
`OPTIONAL_DEPS` list, so the rule fires the moment one is reintroduced; add the
package to the `develop` extra in the same change.

`pylaplace` is a hard dependency and is imported unguarded.

---

## 7. Marker discipline and timeouts

Every test file MUST begin with:

```python
import pytest

pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

Budgets: `unit` -> `timeout(30)` (target < 100 ms per test), `smoke` -> `timeout(60)` (target < 30 s), `integration` -> `timeout(300)`, `slow` -> `timeout(3600)`.

PR CI runs `pytest -m "(unit or smoke) and not skip and not slow and not integration"`; the nightly runs `(integration or slow) and not skip`. Per-function markers are additive, not a replacement. The seed-ensemble statistics carry `integration` and run nightly; a full-system single run with a small system fits `unit` or `smoke` depending on wall time.

---

## 8. Float and numerical comparison

- NEVER `==` for floats; `pytest.approx` or `np.testing.assert_allclose`.
- State the tolerance rationale when non-obvious ("rel=0.2 because the paper quotes a 20 percent scatter against N-body").
- Closed-form algebra pins at `rel=1e-12`; Monte Carlo ensemble statistics at the anchor's stated scatter.

---

## 9. Voice rule for test artifacts

The repo-wide voice rule (zero AI-process disclosure in any public artifact) applies to test code with the same strictness as to source. In scope: test-skip reasons, test-file and test-function docstrings, test names, parametrize ids, log-capture assertions, commit messages on test-touching commits, PR titles and bodies, CI job and step names, inline `src/morrigan/**` comments, shipped log strings. Out of scope: this file, `morrigan-code-review.md`, `copilot-instructions.md`.

Banned inside in-scope artifacts: "audit", "review pass", "adversarial review", AI-roadmap labels, `claude-config/...` paths, "Generated with Claude", AI-tool names, em-dashes, en-dashes (except bibliographic page ranges), process meta-commentary. Write the OUTCOME, never the PROCESS. First-person voice.

---

## 10. Fixture and parameter conventions

- SI units internally (kg, m, s); the config boundary uses Earth masses, AU, and degrees, converted once in `driver.py`. Tests state which side of the boundary their inputs live on.
- Use `@pytest.mark.parametrize` for multiple physical regimes; ids read as scenarios, not tuples.
- Seeds per Section 4. Hypothesis tests use `@settings(derandomize=True)`.
- Use `tmp_path` for temporary files (the TOML settings path).

---

## 11. Documentation per test

- **File-level docstring**: names the source under test and the invariants exercised. Required.
- **Function-level docstring**: physical scenario or contract clause, plus the seed when stochastic. Required (lint-enforced).
- **Inline comments**: why this input range ("masses spanning a decade so the `a^(3/2)` and Hill `m^(1/3)` scalings are resolved above tolerance").

---

## 12. Naming

- Test names describe behavior: `test_kepler_period_scales_as_a_to_three_halves`, NOT `test_kepler_period`.
- Test file names mirror source 1:1. Documented exception: `tests/test_ensemble_statistics.py`, the cross-cutting seed-ensemble statistics exercising the whole model through `run_system` (integration tier).

---

## 13. Independent review trigger

A PR that adds or substantially modifies > 50 lines of test code across all its commits triggers an independent review pass before merge. PR-level denominator: `git diff origin/main...HEAD -- 'tests/**'`. The reviewer cites the anti-happy-path rule, the discrimination-guard requirement (including the argument-order rule from Section 2), the physics-invariant tier, and the seeding discipline from Section 4.

---

## 14. Tooling

- `bash tools/validate_test_structure.sh` -- marker presence and file naming.
- `python tools/check_test_quality.py --check` -- AST scan, blocking on PR.
- `python tools/check_test_quality.py --baseline` -- regenerate the floor after a deliberate sweep.
- `python tools/check_test_quality.py --reference-pinned-status` -- punch list.
- `python tools/update_coverage_threshold.py` -- one-way ratchet, capped at 90.
- `ruff check src/ tests/ tools/` and `ruff format src/ tests/ tools/`.

---

## 15. Failure modes to recognize on review

Real patterns from this codebase's history. The lint script catches some; reviewers catch the rest.

| Pattern | Example | Why it slipped | Fix |
|---|---|---|---|
| **Wrong mass denominator in a fitted ratio** | The loss law evaluated `M_i / M_t` where the fit prescribes `M_i / (M_i + M_t)`; the error is invisible at small mass ratios and ~11-20 percent at the ratios this model produces | The docstring wrote the correct formula while the code evaluated the wrong one; no asymmetric-pair pin existed | Pin an asymmetric pair with the wrong-denominator value in the guard |
| **`exp` vs `10**` in a fitted timescale** | A crossing-time fit published as a base-10 power law evaluated with `np.exp` | Both give "a big number"; no absolute pin at a tabulated point | Pin one tabulated point absolutely; guard with the wrong-base value |
| **Wrong pair index in a mass product** | A stirring term used `Mp[0]*Mp[2]` where the interacting pair was `[0, 1]` | Symmetric test systems (equal masses) hide index errors | Unequal masses across the system in at least one test |
| **Unseeded draw in a test** | A test asserting on stochastic output without seeding | Passed locally by luck of import order | Seed explicitly per Section 4; state the seed in the docstring |
| **Difference-based argument-order guard** | `abs(f(a) - f(b)) > delta` claimed to pin the argument mapping | A permuted dispatch permutes the values and preserves the difference | Absolute pins on both sides plus a direction assertion |
| **Silent skip in helper** | `if actual is None: continue` masking broken introspection | Helper hides a real failure as a no-op | Hard assertion with a message |

When you spot a new variant, add it here.

---

## 16. Sister rules (cross-link)

- `.github/copilot-instructions.md` "Testing Standards" -- the high-level summary.
- `.github/.claude/rules/morrigan-code-review.md` -- the review-pass gate and the domain-aware source-review checks (orbital bounds, unit boundaries, the PROTEUS record-schema contract, global-RNG mutation).

Any change to the rule set: update both files in the same commit and call out the cross-reference in the commit body.
