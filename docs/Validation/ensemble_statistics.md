# Validation: ensemble statistics

Tests: `tests/test_ensemble_statistics.py` (integration tier, nightly).

| Test | Pin | Authority |
| --- | --- | --- |
| `test_eq37_coefficient_matches_the_published_value` | e_ij^2 / e_esc^2 = 0.0688 at b_H = 10, 1 au, 3 g cm-3, 1 Msun, built from the model's own primitives; mass-independent across a hundredfold sweep | Kimura et al. (2025), eq. 37 (published coefficient 0.07) |
| `test_s0_like_ensembles_reach_a_widened_hill_stable_endstate` | Twenty seeds of the equal-mass S0-like setup (15 embryos, 2.43 Earth masses, from 0.1 au at spacing 10) all accrete, widen their mean Hill-scaled separations beyond 10, and settle at moderate Hill-scaled eccentricities; a bounded-window check, not a paper-value pin | Setup from Kimura et al. (2025), Table 1 (model S0); statistics computed per eqs. 33-34 |

The eq. 37 pin is the closed-form anchor: the published 0.07 coefficient is recovered at its stated precision from the eq. 36 crossing eccentricity over the mutual escape eccentricity, and the analytic mass cancellation and quadratic spacing dependence are asserted exactly. The ensemble bounds are physically derived windows, not figure reads: the separation floor is the initial spacing the endstates must exceed to be Hill-stable, and the eccentricity band brackets the excitation the escape-eccentricity scale of the region sets.
