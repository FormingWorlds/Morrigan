# Morrigan Code Review Criteria

When reviewing Morrigan code (either your own or via code-reviewer agents), apply these domain-specific checks in addition to standard code quality review.

> **Discovery note.** Morrigan keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `morrigan-tests.md` explicitly. **Before opening any review pass, read both this file and `morrigan-tests.md`.**

## Physics plausibility

- Masses, radii, densities, and semi-major axes must be strictly positive everywhere. Flag any code path where a merger, a scattering, or the sorting step could produce a non-positive value.
- Eccentricities must stay in `[0, 1)`. Flag any orbital update that can push `e` to or past unity without an explicit hyperbolic-orbit treatment (there is none in this model; such a state is a bug).
- The impact parameter is `b = sin(beta)` and lives in `[0, 1]`.
- The collision speed at contact is `sqrt(v_inf^2 + v_esc^2)` and therefore never below the mutual escape speed. Flag any speed computed another way.
- Merged masses are exact sums. Flag anything that reduces a body's mass in a merger; the model sheds nothing and the coupled framework validates that chain to machine precision.
- The planet count is non-increasing over a run; a merger removes exactly one body.
- Timeline records are strictly ordered in time, and successive records for one body chain their masses.

## Unit convention boundaries

Morrigan is SI internally with a thin conversion layer at the edges:

- **Internal state**: kg, m, s, radians where angles appear.
- **Config boundary** (`driver.py`): Earth masses, AU, degrees, Gyr. Converted once on read; flag any second conversion site.
- **Record boundary** (`run_system` output): SI everywhere except `time` in years, using the model's own 365-day year (`gyr2sec / 1e9`). Flag any record field built from a different year convention.
- **Constants**: single source `src/morrigan/constants.py`. A physical constant introduced as a literal in a function body is a red flag.

## Global random state

The Monte Carlo layer draws from numpy's global random state, seeded once at the `run_system` entry:

- Library code must NOT call `np.random.seed` anywhere except the documented entry-point seeding. A hidden reseed inside a helper silently correlates ensemble members and destroys the seed-to-timeline determinism contract.
- New stochastic code paths must draw through the global state (consistent with the rest of the model) and must keep the same-seed-same-timeline determinism test green.
- Flag any use of `random` (the stdlib module) alongside `np.random`; a second, unseeded RNG breaks determinism invisibly.

## The PROTEUS record-schema contract

`run_system` returns the per-impact records PROTEUS's accretion coupling consumes. The schema (field names, units, and conventions) is documented in `.github/copilot-instructions.md` and is an interface, not an implementation detail:

- Any field rename, unit change, or convention change (contact speed, perfect-merger mass, bulk radii) is a breaking change for PROTEUS and must be called out in the PR description.
- `v_esc` in the record is recomputed from the pair's masses and radii with the mutual formula; keep it consistent with the speed convention.
- Densities in the record are recovered from mass and radius so each record is exactly self-consistent; flag any change that writes an independent density.
- The model carries no atmosphere, so `M_merged_after` is always the plain sum of the two bodies and nothing is shed between impacts. Impact erosion belongs to the consumer: in a coupled run PROTEUS applies `zephyrus.collision.mass_loss` to its own atmosphere, from the geometry the record carries. A change that reintroduces atmosphere tracking here would break that separation and the exact mass chain PROTEUS validates against.

## Fitted-law fidelity

Several formulas are fits transcribed from Kimura et al. (2025). For any transcription change:

- Verify the base of every exponential against the paper (`10**` vs `exp` is a historical bug class here).
- Verify every ratio's denominator against the paper (`M_i / M_tot` vs `M_i / M_t` is a historical bug class here).
- Verify the fit's validity domain is stated in the docstring, and that any clamp needed to keep the output physical is present.
- A transcription PR must cite the equation number in the paper.

## Config mutability

The parsed settings dictionary carries user input and is not mutated after the run starts. Flag any code that writes back into it mid-run; use local variables.

## Test marker discipline

Every test file must begin with a module-level `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]` (unit/30 s, smoke/60 s, integration/300 s, slow/3600 s). Per-function markers are additive; CI runs the fast tiers on PR and the slow tiers nightly, and any file missing the tier marker ships untested.

## Test quality (cross-reference)

Test-content rules (anti-happy-path, discriminating-value guards including the argument-order rule, physics-invariant tiering, certification markers, seeding discipline, mocking discipline) live in [`morrigan-tests.md`](morrigan-tests.md). When reviewing tests, apply both files: this one for marker discipline and the source-side checks, the deep-dive for the content contract.

## Sister rules (cross-link)

- [`.github/copilot-instructions.md`](../../copilot-instructions.md) "Testing Standards" -- high-level rules visible to all readers. Repo-root `CLAUDE.md` is a symlink to this file.
- [`morrigan-tests.md`](morrigan-tests.md) -- test quality deep-dive.
