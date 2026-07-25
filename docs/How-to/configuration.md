# Configuration

Morrigan reads one TOML settings file, passed with `-c`. `initialise.toml` in the repository root is a worked example. This page covers what each block does and how the settings interact; [Settings and inputs](../Reference/parameters.md) is the field-by-field reference.

## The three blocks

```toml
['run_simulation']
    t              = 0.0            # leave at zero for a fresh run
    t_ref          = 0.0            # leave at zero for a fresh run
    t_event        = 0.0            # leave at zero for a fresh run
    flag_event     = 1              # leave at 1
    a_min          = 0.005
    max_time       = 1
    random_seed    = 1
    save_directory = '1au_system'

['batch']
    ndisk = 20
    nproc = 8

['init_par']
    N          = 10
    Mp         = [0.1, 1.5, 0.9, 0.6, 0.75, 1.2, 1.4, 0.45, 0.5, 1.0]
    e          = 0.01
    inner_edge = 1
```

The `init_par` block is an excerpt; the other two are complete. [Settings and inputs](../Reference/parameters.md) is the field reference, with every field, its meaning and its units; this page is about how the settings interact.

Every key in `run_simulation` is required and read directly, so omitting one fails with a bare `KeyError` naming it. `random_seed` is the only one with a default. The first four hold the clock and the secular bookkeeping and are always zero, zero, zero and 1 for a fresh run; `initialise.toml` in the repository root carries them with the same values.

`Mp` must have `N` entries. `spacing` is optional and defaults to 10 mutual Hill radii.

## Choosing the initial conditions

**Spacing and eccentricity set the instability time**, and it depends steeply on both. A widely spaced, nearly circular system may not go unstable at all inside the evolution time, in which case the run finishes with no mergers. If that happens, tighten the spacing or raise the initial eccentricity before reaching for a longer `max_time`.

There is an upper limit on spacing. The layout condition has a pole where the requested gap approaches the span it is measured across, roughly where `spacing * ((M1+M2)/3M*)^(1/3)` reaches 2, and beyond it the layout is refused with an error naming the pair. It is reachable with giant-planet masses at wide spacings, not with terrestrial ones.

**The mass distribution shapes the outcome more than the total.** Equal-mass systems behave differently from staggered ones: the excitation a scattering imparts is set by the ratio of the two masses, so equal-mass pairs share their kick evenly and unequal pairs throw the lighter body much harder. An equal-mass system is the cleaner comparison against published statistics; a staggered one is closer to a real disk.

**`impact_angle` is a single value for the whole run.** The model does not draw a distribution of impact geometries, so every collision in a run happens at the same angle. Its sine is the impact parameter, so 90° is a grazing hit and 0° is head-on.

## Seeds and reproducibility

Each system's seed is mixed from `random_seed` and the system's index within the batch. Two consequences:

- A settings file reproduces its systems exactly, however many times it is run and whatever `nproc` is set to.
- Two batches with neighbouring `random_seed` values are independent. Estimating the scatter between ensembles by rerunning with a few different seeds is therefore a valid thing to do.

## Batch size and parallelism

`ndisk` is how many systems to evolve, `nproc` how many worker processes to spread them over. A single system takes seconds, so `ndisk` is usually set by how many realisations the statistics need rather than by runtime. `nproc` above the core count gains nothing.

## Output location

Results are written under `save_directory`. A relative path is taken from the directory the command runs in, not from where the settings file lives, so the run prints the full path it resolved. Keep the settings file alongside the results it produced; the repository ignores `*.toml` other than `initialise.toml`, so a settings file worth keeping has to be added to git explicitly.

## Driving PROTEUS instead

None of this applies when Morrigan runs inside PROTEUS: the settings file is not read, and the same quantities come from the `[accretion.morrigan]` block of the PROTEUS configuration. See [Coupling to PROTEUS](proteus_coupling.md).
