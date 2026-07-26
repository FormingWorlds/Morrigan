# Tutorial: run PROTEUS with Morrigan

This is one complete PROTEUS configuration that grows a planet by giant impacts,
with every other module set to `dummy` so it runs in seconds and needs no
external solver. It is the shortest way to see the coupling work before putting
Morrigan in front of a real interior and atmosphere.

!!! info "What this assumes"
    The accretion module is part of PROTEUS from
    [pull request #789](https://github.com/FormingWorlds/PROTEUS/pull/789)
    onwards. Use a PROTEUS that has it: either `main` once that work is merged,
    or a checkout of its branch. Morrigan itself is an optional dependency:

    ```bash
    pip install "fwl-proteus[morrigan]"
    ```

## The configuration

Save this as `morrigan_tutorial.toml`:

```toml
# PROTEUS: all-dummy configuration
# Every module set to "dummy" for fast test runs without external solvers
# (no AGNI, SPIDER, VULCAN, MORS, Zalmoxis, CALLIOPE, or SOCRATES needed).
# Use this for wiring tests, CI, and rapid development iteration.
# Config file format version
config_version = "3.0"

    [params.out]
        path = "auto"

    [params.dt]
        initial  = 1e2           # initial step [yr]
        minimum  = 1e2           # min step [yr]
        maximum  = 3e7           # max step [yr]

        [params.stop.time]
            maximum = 1e9        # stop time [yr]
        [params.stop.solid]
            enabled = true
        [params.stop.escape]
            enabled = false
        [params.stop.radeqm]
            enabled = false

# Star: fixed solar luminosity (no evolution).
# Morrigan reads the host mass from here rather than repeating it.
[star]
    module = "dummy"
    mass   = 1.0                 # [M_sun]
    [star.dummy]
        radius = 1.0                 # [R_sun]

# Orbit: 0.5 AU with weak tidal heating
[orbit]
    module        = "dummy"
    semimajoraxis = 0.5          # [AU]
    eccentricity  = 0.1
    [orbit.dummy]
        H_tide = 1e-9            # tidal power density [W/kg]
        Imk2   = -0.01           # Im(k2) love number

# Planet: 1 Earth mass
[planet]
    mass_tot      = 1.0              # [M_earth]
    # adiabatic_from_cmb keeps the all-dummy path free of the silicate
    # liquidus lookup, so this config runs without any external solver.
    temperature_mode = "adiabatic_from_cmb"
    tsurf_init       = 3000.0        # [K], the state each impact re-melts to
    volatile_mode = "elements"
    [planet.elements]
        O_mode   = "ic_chemistry"
        O_budget = 0.0
        H_mode   = "ppmw"
        H_budget = 1e4
        C_mode   = "ppmw"
        C_budget = 1e3
        N_mode   = "ppmw"
        N_budget = 500
        S_mode   = "ppmw"
        S_budget = 500

# Interior structure: scaling laws (no external solver)
[interior_struct]
    module         = "dummy"
    core_frac      = 0.55        # radius fraction
    core_frac_mode = "radius"

# Interior energetics: parameterized cooling
[interior_energetics]
    module = "dummy"
    [interior_energetics.dummy]
        mantle_tliq = 2700.0
        mantle_tsol = 1700.0
        mantle_cp   = 5000.0     # higher heat capacity = slower cooling

# Outgassing: parameterized partitioning (no thermodynamic solver)
[outgas]
    fO2_shift_IW = 2
    module       = "dummy"

# Atmosphere: grey opacity model
[atmos_clim]
    module     = "dummy"
    surf_state = "fixed"
    rayleigh   = false
    albedo_pl  = 0.1
    [atmos_clim.dummy]
        gamma = 0.5              # T_rad = T_surf * (1 - gamma); 0 = transparent, 1 = zero OLR

# Escape: fixed bulk rate
[escape]
    module = "dummy"
    [escape.dummy]
        rate = 0.0               # no escape, so the run reaches solidification

[atmos_chem]
    module = "dummy"
    when   = "offline"

# Accretion: grow the planet by giant impacts from a Morrigan system
[accretion]
    module             = "morrigan"
    time_offset        = 0.0         # Morrigan's clock starts at disk dispersal
    impactor_volatiles = "match_planet"
    atmloss_module     = "none"

    [accretion.morrigan]
        seed              = 7
        num_planets       = 6
        masses            = [0.4, 0.8, 0.5, 1.0, 0.6, 0.7]   # [M_earth]
        eccentricity_init = 0.01
        inner_edge        = 0.05     # [AU]
        spacing           = 10.0     # mutual Hill radii
        density           = 5500.0   # [kg m-3]
        impact_angle      = 45.0     # [deg]
        evolution_time    = 0.5      # [Gyr]
        inner_cutoff      = 0.005    # [AU]
        selector          = "mass"
```

## Run it

```bash
proteus start -c morrigan_tutorial.toml --offline
```

## What happens

Morrigan runs once, before the loop, and hands PROTEUS the impact history of one
survivor. The start-up lines report the system it evolved and the schedule it
kept:

```
[ INFO  ] Accretion module  morrigan
[ INFO  ] Preparing accretion model
[ INFO  ] Running giant-impact model for 6 embryos
[ INFO  ] Following body 1 (selector 'mass'): 0.800 -> 1.700 M_earth, 0.0556 -> 0.0558 AU
[ INFO  ] Body 1 experienced 2 impacts
[ INFO  ] Scheduled 2 impact(s)
[ INFO  ]     first at 2.6203e+03 yr, last at 3.1016e+03 yr
```

Each impact is then applied as the loop reaches its time:

```
[ INFO  ] Giant impact at t = 2.6203e+03 yr: target 1 struck by 2, adding 0.5000 M_earth
[ INFO  ]     planet is now 1.4964 M_earth at 0.52073 AU, e = 0.0858
[ INFO  ] Giant impact at t = 3.1016e+03 yr: target 1 struck by 0, adding 0.4000 M_earth
[ INFO  ]     planet is now 1.8934 M_earth at 0.50202 AU, e = 0.0502
[ INFO  ] ===> Planet solidified!! <===
```

The run finishes at 2.94 &times; 10<sup>4</sup> yr with the mantle solid. Two
things in those lines are worth understanding.

**The planet grew, and its orbit moved.** It started at 1.0 M&#8853; on a 0.5 AU
orbit with an eccentricity of 0.1, and ended at 0.502 AU with an eccentricity of
0.0502. The dynamical system sits at 0.056 AU, nowhere near 0.5 AU, because only
the *change* each impact made is transferred: the semi-major axis takes the
ratio and the eccentricity the difference. The planet keeps the orbit you
configured and inherits the kick.

**The mass in the log is not the whole planet.** The `planet is now ... M_earth`
line reports the interior anchor, the rock. The helpfile's `M_planet` also
carries the volatile budgets, so it ends at 1.961 M&#8853; while
`M_accreted_rock`, the cumulative rock the impacts delivered, reads 0.893
M&#8853;.

The run also writes the schedule it resolved to `impact_timeline.csv` in the
output directory. Resuming replays that file, and
[`accretion.module = "timeline"`](../How-to/timeline_module.md) replays it as a
run of its own.

## Things to change

- **`accretion.morrigan.inner_edge`** is what makes this tutorial quick. A
  compact system at 0.05 AU goes unstable within a few thousand years, so the
  impacts land while the magma ocean is still there. Move it out to 1 AU and the
  instability time grows past the point where this all-dummy planet has already
  solidified, and the run ends with the impacts still ahead of it. That is not a
  failure of the coupling; it is the schedule and the interior evolving on
  different clocks, and it is the first thing to check when a run reports
  scheduled impacts but never applies one.
- **`accretion.impactor_volatiles`** is `match_planet` here, so each impactor
  carries the planet's own formation abundances. `dry` adds rock only.
- **`accretion.atmloss_module`** is off. Set it to `zephyrus` for the
  Kegerreis et al. (2020) erosion law, or `constant` with `atmloss_frac`.
- **`interior_energetics.module`** is `dummy` for speed. `aragog` is the
  production choice and the only one that books the re-melt heat into the energy
  budget; `spider` is refused with accretion enabled, because it has no
  validated re-melt path.

## Next

- [Coupling to PROTEUS (how-to)](../How-to/proteus_coupling.md) for the full
  recipe, the selector choice, and the pitfalls.
- [Coupling to PROTEUS (explanation)](../Explanations/proteus_coupling.md) for
  what each impact does and why the coupling is built this way.
- [Run a system standalone](standalone.md) to see the dynamics on their own.
