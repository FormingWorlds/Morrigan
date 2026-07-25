# Validation: merger bookkeeping

Source: `src/morrigan/merge_embryo.py`. Tests: `tests/test_merge_embryo.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_collision_velocity_floors_at_the_mutual_escape_speed` | v_c equals the mutual escape speed exactly for circular orbits; the eccentric case pinned at 11945.7 m/s | The sqrt(v_inf^2 + v_esc^2) contact-speed convention (analytical limit) |
| `test_merger_closes_its_mass_and_orbit_books` | The merged orbit lands at (M1+M2)/(M1/a1 + M2/a2), 1.05882 au for the reference pair | Kimura et al. (2025), eq. 16 |

Guards: the wrong linear combination v_inf + v_esc, and the arithmetic-mean and mass-weighted-mean wrong orbits. Exact mass closure on the plain sum, mass-based survivor selection, and the eq. 17 eccentricity composition bounds are asserted alongside.
