# Validation anchors

Every physics source file in Morrigan carries at least one test pinned against a published benchmark, an analytical limit, or a cross-implementation cross-check, tagged `@pytest.mark.reference_pinned`. The pages in this section are the per-source inventory: what is pinned, against what authority, and in which test. The companion marker `@pytest.mark.physics_invariant` tags every test that asserts a physical invariant (conservation, positivity or boundedness, monotonicity or symmetry, or a discriminating pinned value); the reference-pinned tests are the subset whose authority is external to the implementation.

| Source file | Anchor | Authority |
| --- | --- | --- |
| [helper_functions.py](helper_functions.md) | Kepler's third law, the Earth Hill radius, the Earth mean radius | Analytical limits and published values |
| [interaction_timescales.py](interaction_timescales.md) | Timescale closed forms; the analytic stability boundary | Kimura et al. (2025) eqs. 6-11, 23, 28 |
| [crossing_pair.py](crossing_pair.md) | Perihelion-gap pair selection | Kimura et al. (2025) eq. 4 |
| [secular_solution.py](secular_solution.md) | Epoch recovery of the eigen-decomposition | Laplace-Lagrange theory (Murray & Dermott 1999) |
| [merge_embryo.py](merge_embryo.md) | The contact-speed floor; the merged orbit | Analytical limit; Kimura et al. (2025) eq. 16 |
| [orbit_cross_K25.py](orbit_cross_K25.md) | Scattering conservation of the mass-weighted orbit sum | Kimura et al. (2025) eqs. 18-20 |
| [driver.py](driver.md) | Record self-consistency against the mutual escape speed | Analytical limit |
| [Ensemble statistics](ensemble_statistics.md) | The eq. 37 coefficient; Hill-scaled endstates of the S0-like ensemble | Kimura et al. (2025) eqs. 33-37, Table 1 |

The punch list of sources still lacking a pin is reported by `python tools/check_test_quality.py --reference-pinned-status`; it is empty at the current state of the suite.
