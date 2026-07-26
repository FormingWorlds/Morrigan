# Tutorial: run a system standalone

This runs one system of six embryos from a settings file, in under a second, and
reads the result. It needs Morrigan installed and nothing else; see
[installation](../How-to/installation.md) if you have not done that yet.

## The settings file

Save this as `tutorial.toml`:

```toml
['run_simulation']
# Bookkeeping the driver expects; leave these four as they are.
t = 0.0
t_ref = 0.0
t_event = 0.0
flag_event = 1

a_min = 0.005            # perihelion inside which a body is lost to the star [AU]
max_time = 0.5           # evolution time [Gyr]
save_directory = 'tutorial_run'
random_seed = 7          # fixing this reproduces the run exactly

['batch']
ndisk = 1                # one system, so the output is easy to read
nproc = 1

['init_par']
N = 6                    # six embryos
e = 0.01                 # shared initial eccentricity
Mp = [0.4, 0.8, 0.5, 1.0, 0.6, 0.7]   # embryo masses [M_earth]
impact_angle = 45        # [deg]

Ms = 1                   # stellar mass [M_sun]
rho_p = 5500             # bulk density [kg m-3]
inner_edge = 0.3         # innermost embryo [AU]
```

The section header on the first line is required: the driver reads its settings
from `['run_simulation']`, and a file without it stops with a `KeyError`. Every
field is described in [settings and inputs](../Reference/parameters.md).

## Run it

```console
$ morrigan -c tutorial.toml
Writing results to /path/to/tutorial_run
Ran 1 system in 0.04s
```

Four files appear under `tutorial_run`:

| Path | Contents |
| --- | --- |
| `data/full_systems/full_system_00.csv` | every body's orbit and mass through time |
| `data/mergers/mergers_00.csv` | one row per collision |
| `data/survivors/survivors_00.csv` | what was left at the end |
| `batch_summary.csv` | runtime and surviving-body count |

## What this system does

Six embryos between 0.28 and 0.52 AU go unstable early. Four collisions occur,
all within the first 190 000 yr, and two bodies survive:

| Time \[yr] | Target | Impactor | Target before \[M⊕] | Merged after \[M⊕] |
| --- | --- | --- | --- | --- |
| 60 120 | 1 | 0 | 0.80 | 1.20 |
| 82 458 | 4 | 2 | 0.60 | 1.10 |
| 90 025 | 4 | 3 | 1.10 | 2.10 |
| 187 332 | 4 | 1 | 2.10 | 3.30 |

```
| id |  Mp [kg]    | a_AU   | ecc   |
|  4 | 1.97129e+25 | 0.3335 | 0.074 |
|  5 | 4.18152e+24 | 0.7394 | 0.394 |
```

Body 4 grows from 0.60 to 3.30 M⊕ in three collisions and settles at 0.33 AU on
a nearly circular orbit. Body 5 never collides: it is scattered from 0.52 AU out
to 0.74 AU and left with an eccentricity of 0.39, which is what a system does to
the body that survives without merging.

## The result

<figure markdown>
  ![Semi-major axis and mass of each body against time](../assets/tutorial_standalone_light.png#only-light)
  ![Semi-major axis and mass of each body against time](../assets/tutorial_standalone_dark.png#only-dark)
  <figcaption>
    Each body's semi-major axis (top) and mass (bottom) against time, for the
    settings file above. A track that ends in a dot is a body that was merged
    away; the step changes in the surviving track are its collisions.
  </figcaption>
</figure>

Two things are worth reading off it. The whole giant-impact phase is over within
the first fraction of a percent of the evolution time, which is why the tracks
are flat across most of a logarithmic axis: once the system is Hill-stable
nothing further happens, and running longer changes nothing. And the mass that
body 4 gains is exactly the mass the other tracks lose, because a collision in
this model is a perfect merger.

The figure above is drawn by `tools/plot_tutorial.py` in the repository. The
`plot.py` script at the repository root plots a whole batch instead, with
per-system tracks, orbits and ensemble statistics:

```console
$ pip install -e ".[plot]"
$ python plot.py -c tutorial.toml
```

## What to change next

- `random_seed` samples a different realisation of the same initial conditions.
  One seed gives one history, bit for bit; the spread across seeds is the
  meaningful quantity, not any single run.
- `N` and `Mp` set how much material the system has and how it is divided.
- `inner_edge` moves the whole system in or out, which changes the collision
  speeds and the instability time.
- `ndisk` in `['batch']` runs many systems on a process pool for the statistics.

[Choosing the initial conditions](../How-to/configuration.md#choosing-the-initial-conditions)
covers how spacing and eccentricity set the instability time, which is what
decides whether a run has any collisions at all.

## Next

- [Run Morrigan inside PROTEUS](proteus_minimal.md) with a complete configuration file.
- [Model overview](../Explanations/model.md) for the physics behind these tracks.
- [Limitations](../Explanations/limitations.md) for what this model does not represent.
