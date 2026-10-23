# Validation: impact melt

Source: `src/morrigan/impact_melt.py`. Tests: `tests/test_impact_melt.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_port_equals_the_reference_code` | f_melt, f_melt_raw, its 1-sigma, dIE and h equal the reference to 1e-12 relative on 120 cases: angles 0, 30, 45, 60 and 90°, mass ratios 0.03, 0.1, 0.3 and 0.5, speeds 1 and 1.5 v_esc, total masses 2, 9.31 and 40 Mars masses; the 42 cases the reference fits as hit-and-run are compared with its merging branch forced and flagged above v_cr | Cross-implementation check against MeltScalingLaw at commit 2eca30f2 (Nakajima et al. 2021), fed its own radii and Morrigan's G |
| `test_issue_table_at_one_earth_mass` | The 9 cells of the melt table at 1 M_E and v_esc to two decimals; the 90°, 0.3 cell's 1-sigma is 0.38 in the reference code | The same reference |

`tools/generate_nakajima_golden.py` writes `tests/data/nakajima2021_golden.json` from a checkout of the reference, patched in memory only (`np.float`, colour maps, G, an optional forced merging branch, and a return of the bulk values before the clip). The reference code's sixth-order Legendre term is twice the standard P6 (1/8 in place of 1/16), and the port evaluates that basis, the one its coefficients come with; with the standard P6 the internal-energy fit differs by up to 2.7 % on these cases.

Guards: the latent heat in place of E_M, the 0.70 mantle fraction left out, and dIE without its angle fit each land far from the pinned value. Monotonic decline from head-on to grazing over 181 angles, the bounds of every flag, the merged radius at the target density, the 1-sigma at each exact tie between simulated angles, and the refusal of every argument outside its domain are asserted alongside; the Hypothesis property tests are in `tests/test_impact_melt_properties.py`.
