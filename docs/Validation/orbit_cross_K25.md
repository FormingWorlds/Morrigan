# Validation: crossing outcomes

Source: `src/morrigan/orbit_cross_K25.py`. Tests: `tests/test_orbit_cross_K25.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_scattering_conserves_the_mass_weighted_orbit_sum` | sum(M a) conserved to machine precision through a scattering, with the inner orbit moving in and the outer out | Kimura et al. (2025), eqs. 18-20 |

The conservation is an identity of the published orbit shifts, so any transcription error in the mass weighting breaks it at the percent level. The merger path is pinned to its record schema and mass books, the ejection path to removing exactly the body excited past e = 1 with the survivor bound inside its original orbit, and the stable-pair path to complete inaction, including not consuming a random draw.
