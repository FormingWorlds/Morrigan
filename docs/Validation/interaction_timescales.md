# Validation: interaction timescales

Source: `src/morrigan/interaction_timescales.py`. Tests: `tests/test_interaction_timescales.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_tau_vis_pins_its_closed_form` | 1.74771e11 s for two Earths at 0.95/1.05 au, e = 0.05; exact rep_e^4 scaling | Kimura et al. (2025), eqs. 8-10, 23 |
| `test_tau_vis_pins_its_closed_form` (unequal-mass case) | 6.41894e9 s for 1 and 10 Earth masses at e = 0.001 | Kimura et al. (2025), eq. 6 |
| `test_tau_col_pins_and_returns_a_time_not_a_rate` | 1.97470e12 s for the equal-mass pair; 2.43544e11 s for 1 and 10 Earth masses at e = 0.001 | Kimura et al. (2025), eqs. 6, 11 |
| `test_interaction_wrapper_pins_the_analytic_stability_boundary` | Stability flips inside the bracket \[4.95, 5.10\] around the analytic root x_crit = 5.0332 for e = 2h | Kimura et al. (2025), eq. 28 |

Eq. 6 is resolved only by the unequal-mass pins. Eq. 23 takes the larger of `sum(e_cross)` and `sum(e)`, and for the equal-mass reference pair those agree to within one bit, so the crossing eccentricities do not set `rep_e` there. Both timescales carry the shared denominator, so each is pinned separately: the swapped-denominator variant lands 4.8 % away for the viscous timescale and 2.4 % away for the collision timescale.

Guards: the factor-3 variant and dropped-2pi variant of the viscous timescale, the forgotten rate inversion of the collision timescale (twenty-five orders of magnitude), and the wrong-constant variants of the stability criterion, whose roots (4.62 and 5.77) fall outside the bracket. The Petit et al. (2020) crossing time is ordered across 4, 6, and 8 mutual Hill radii, with its overlap (0.0) and stable (1e20) sentinels asserted exactly.
