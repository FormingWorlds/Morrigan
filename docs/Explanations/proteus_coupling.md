# Coupling to PROTEUS

PROTEUS consumes Morrigan through the in-memory entry point `morrigan.run_system`, which evolves one system and returns its survivors and the full impact history of each, without writing any files. The per-impact record schema is the interface between the two codes.

## The call

```python
from morrigan import run_system

out = run_system(
    seed=42,
    masses=[...],            # embryo masses [kg]
    eccentricity=0.05,
    inner_edge=0.05 * 1.5e11,  # [m]
    spacing=10.0,              # [mutual Hill radii]
    density=5500.0,            # [kg m-3]
    impact_angle=45.0,         # [deg]; its sine is the impact parameter
    evolution_time=1.0,        # [Gyr]
    inner_cutoff=0.005 * 1.5e11,  # [m]
    stellar_mass=1.0,          # [Msun]
)
```

Inputs are SI apart from the stellar mass (solar masses) and the evolution time (Gyr). The returned quantities are SI with times in years, using the model's own 365-day year.

## The impact record schema

`out['impacts']` is a dictionary keyed by surviving-body id; each value is that body's impacts in time order, and every surviving body is present, with an empty list if it never merged. Each record carries:

`time [yr]`, `M_target_before [kg]`, `M_impactor [kg]`, `M_merged_after [kg]`, `v_impact [m/s]`, `v_esc [m/s]`, `impact_parameter`, `R_target_before [m]`, `R_impactor [m]`, `rho_target [kg/m3]`, `rho_impactor [kg/m3]`, `a_before [m]`, `a_after [m]`, `e_after`, `id_target`, `id_impactor`.

Conventions the schema guarantees:

- `v_impact` is the speed at first contact, $\sqrt{v_\infty^2 + v_\mathrm{esc}^2}$, the convention the Kegerreis et al. (2020)[^cite-kegerreis2020] erosion law expects, and `v_esc` is the mutual escape speed recomputed from the recorded masses and radii, so the two are exactly consistent.
- `M_merged_after` is the perfect-merger sum of the two masses. The model sheds nothing, so a consumer applying its own erosion law starts from the unstripped mass.
- Bodies carry no atmosphere at all, in the record or in the model; masses and radii are bulk values, and the densities are recovered from mass and radius, so each record is exactly self-consistent.
- Successive records for one body chain: each `M_target_before` equals the previous record's `M_merged_after` to machine precision, since nothing removes mass between impacts.

`out['survivors']` lists each surviving body's id together with its initial and final mass and semi-major axis, so the consumer can select a planet and replay its growth.

Any change to the record fields, units, or conventions is a breaking interface change for PROTEUS and must be flagged in the pull request that makes it.

## Random-state hygiene

`run_system` seeds numpy's global random state for its own draws and restores the caller's state afterwards, so an embedding program's random sequence is unaffected by running a system in-process.

## References

[^cite-kegerreis2020]: Kegerreis, J.A., Eke, V.R., Catling, D.C., Massey, R.J., Teodoro, L.F.A. & Zahnle, K.J., *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020.
