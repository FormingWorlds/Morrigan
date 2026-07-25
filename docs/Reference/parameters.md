# Settings and inputs

## The TOML settings file

The command line reads one TOML file (`morrigan -c <file>`). Units at this boundary are Earth masses, au, degrees, and Gyr; everything is converted once on read and the model runs in SI.

### `[run_simulation]`

| Setting | Meaning | Units |
| --- | --- | --- |
| `t`, `t_ref`, `t_event` | Initial clock, secular reference epoch, and first event time; all zero for a fresh run | s |
| `flag_event` | 1 to compute the secular solution and first crossing at start; leave at 1 | - |
| `a_min` | Perihelion inside which a body has fallen into the star | au |
| `max_time` | Length of the dynamical evolution | Gyr |
| `random_seed` | Base seed; system i runs with `random_seed + i` | - |
| `save_directory` | Where result tables are written; relative paths resolve from the working directory | - |

### `[init_par]`

| Setting | Meaning | Units |
| --- | --- | --- |
| `N` | Number of embryos | - |
| `e` | Initial eccentricity shared by all embryos | - |
| `impact_angle` | Impact angle for every collision; its sine is the impact parameter | deg |
| `Mp` | Embryo masses, one per embryo | Earth masses |
| `Ms` | Stellar mass | solar masses |
| `rho_p` | Bulk density shared by all embryos | kg m-3 |
| `inner_edge` | Orbit of the innermost embryo | au |
| `spacing` | Spacing between adjacent embryos, in Hill radii evaluated at the inner orbit of each pair (see [Model overview](../Explanations/model.md)) | Hill radii |

### `[batch]` (optional)

| Setting | Meaning |
| --- | --- |
| `ndisk` | Number of systems to run (default 1) |
| `nproc` | Worker processes for a batch (default: cores minus one) |

## The in-memory entry point

`morrigan.run_system(seed, masses, eccentricity, inner_edge, spacing, density, impact_angle, evolution_time, inner_cutoff, stellar_mass)` takes SI inputs (kg, m) apart from the stellar mass in solar masses and the evolution time in Gyr, and returns survivors and per-impact records in SI with times in years. The record schema is documented in [Coupling to PROTEUS](../Explanations/proteus_coupling.md).
