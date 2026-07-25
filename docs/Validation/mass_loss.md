# Validation: impact atmosphere loss

Source: `src/morrigan/mass_loss.py`. Tests: `tests/test_mass_loss.py`, `tests/test_mass_loss_properties.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_pins_the_kegerreis_scaling_closed_form` | X = 0.147467 for a 0.1 mass-ratio impactor at 1.5 v_esc, equal densities, b = 0.7 | Kegerreis et al. (2020), eq. 1 with the volume-fraction geometry factor |
| `test_equal_density_agreement_with_zephyrus` | b-dependence ratios agree with `zephyrus.collision.mass_loss` to 1e-9; absolute values to 1.5e-3 | Cross-implementation cross-check |

Guards: the M_i/M_t denominator slip (0.152106) and the 0.5 outer exponent slip (0.206922) both sit more than one hundred tolerances from the pin. The unclamped head-on value 1.25801 documents that clamping is the call site's responsibility, and b = 1 is asserted as an exact algebraic zero. The residual 1.5e-3 band in the absolute cross-check is the (G_morrigan/G_zephyrus)^0.65 factor between the two packages' gravitational constants.
