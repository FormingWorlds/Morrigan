# Validation: helper functions

Source: `src/morrigan/helper_functions.py`. Tests: `tests/test_helper_functions.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_kepler_period_pins_the_analytic_third_law` | P(1 au, 1 Msun) = 3.16894e7 s with the model constants; P(4a)/P(a) = 8 exactly | Kepler's third law (analytical) |
| `test_hill_sphere_reproduces_the_earth_value` | Earth's Hill radius 1.5005e9 m (0.0100 au); 8-fold mass gives exactly 2-fold radius | Published Earth value; the (M/3Ms)^(1/3) closed form |
| `test_planet_radius_and_density_invert_each_other` | R(M_earth, 5513 kg m-3) = 6.3717e6 m; exact round trips at both mass extremes | Published Earth mean radius and density |

Discrimination guards pin the formulas, not just the values: dropping the factor 3 in the Hill denominator (44 percent high), a linear-in-a period (factor 8 vs 4), and the 4/3-vs-3/4 sphere-volume slip (21 percent) all land far outside the tolerances.
