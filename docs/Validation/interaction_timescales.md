# Validation: interaction timescales

Source: `src/morrigan/interaction_timescales.py`. Tests: `tests/test_interaction_timescales.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_tau_vis_pins_its_closed_form` | 1.74771e11 s for two Earths at 0.95/1.05 au, e = 0.05; exact rep_e^4 scaling | Kimura et al. (2025), eqs. 6, 8-10, 23 |
| `test_interaction_wrapper_pins_the_analytic_stability_boundary` | Stability flips inside the bracket [4.95, 5.10] around the analytic root x_crit = 5.0332 for e = 2h | Kimura et al. (2025), eq. 28 |

Guards: the factor-3 variant and dropped-2pi variant of the viscous timescale, the forgotten rate inversion of the collision timescale (twenty-five orders of magnitude), and the wrong-constant variants of the stability criterion, whose roots (4.62 and 5.77) fall outside the bracket. The Petit et al. (2020) crossing time is ordered across 4, 6, and 8 mutual Hill radii, with its overlap (0.0) and stable (1e20) sentinels asserted exactly.
