# Coupling to PROTEUS

This page is the **recipe**: the TOML to write, the choices to make, and the mistakes that cost an afternoon. For why the coupling is built this way, see the [explanation page](../Explanations/proteus_coupling.md).

Morrigan is an optional PROTEUS module. Install it with

```bash
bash tools/get_morrigan.sh
```

from a PROTEUS checkout, which clones the pinned commit and installs it editable.

## Minimal `[accretion]` block

```toml
[accretion]
    module      = "morrigan"
    time_offset = 0.0

    [accretion.morrigan]
        seed              = 1
        num_planets       = 10
        mass_equal        = 0.5      # M_earth, used when masses = []
        eccentricity_init = 0.01
        inner_edge        = 0.1      # AU
        spacing           = 10.0     # mutual Hill radii
        density           = 5500.0   # kg m-3
        impact_angle      = 45.0     # deg
        evolution_time    = 1.0      # Gyr
        inner_cutoff      = 0.005    # AU
        selector          = "match_config"
```

That alone gives a planet that grows by impacts and moves orbit. It delivers no volatiles and strips no atmosphere; both are opt-in below.

The stellar mass is **not** in this block: it comes from `star.mass`. See the [explanation](../Explanations/proteus_coupling.md#from-morrigans-record-to-a-proteus-event) for why.

## Choosing a selector

Morrigan leaves several survivors and PROTEUS follows one. This is the decision that most changes what you get.

```toml
selector = "match_config"    # closest to the configured planet (default)
selector = "mass"            # the most massive survivor
selector = "semimajoraxis"   # nearest a target orbit
selector_value = 1.0         # AU, required for semimajoraxis
selector = "id"              # a specific embryo
selector_value = 3           # embryo index, required for id
```

`match_config` is the sensible default: you have configured a planet mass and orbit, and it finds the survivor that most resembles it. Reach for `mass` when you want whatever the system's dominant body turned out to be, and for `semimajoraxis` when the orbital distance is what you are studying.

`selector_value` is required for `semimajoraxis` and `id` and rejected at config load if missing.

## Delivering volatiles

By default impactors are dry: they add rock and iron only, so the planet's bulk volatile concentration falls by dilution as it grows.

```toml
[accretion]
    impactor_volatiles = "match_planet"
```

`match_planet` gives every impactor the planet's own initial fractional abundances, scaled to the impactor mass, on the assumption that all embryos formed from the same disk material. The composition is frozen at formation, so it does not track the planet's later evolution.

For explicit control:

```toml
[accretion]
    impactor_volatiles = "ppmw"
    impactor_H_ppmw    = 100.0
    impactor_C_ppmw    = 50.0
```

Setting a ppmw budget under any other mode is rejected at config load rather than silently ignored, which otherwise produces a run that looks configured but delivers nothing.

## Stripping atmosphere

Also off by default.

```toml
[accretion]
    atmloss_module = "zephyrus"     # Kegerreis et al. (2020) erosion law
```

evaluates the scaling law from each impact's own collision parameters: speed, mass ratio, density ratio and impact angle. It requires ZEPHYRUS, which PROTEUS installs as standard.

```toml
[accretion]
    atmloss_module = "constant"
    atmloss_frac   = 0.1            # 10 % per impact
```

applies a fixed fraction instead, which is useful for isolating the effect of erosion from its dependence on impact geometry.

One fraction governs both bodies at each impact; the [explanation](../Explanations/proteus_coupling.md#what-proteus-adds-that-morrigan-does-not-do) sets out that convention.

## Worked example: a compact system feeding a growing planet

```toml
[star]
    module = "mors"
    mass   = 1.0                    # Morrigan reads the stellar mass here

[accretion]
    module             = "morrigan"
    time_offset        = 0.0
    impactor_volatiles = "match_planet"
    atmloss_module     = "zephyrus"

    [accretion.morrigan]
        seed              = 7
        num_planets       = 8
        masses            = [0.2, 0.9, 0.4, 1.1, 0.3, 0.7, 1.3, 0.5]
        eccentricity_init = 0.05
        inner_edge        = 0.05
        spacing           = 10.0
        density           = 5500.0
        impact_angle      = 20.0
        evolution_time    = 1.0
        inner_cutoff      = 0.005
        selector          = "mass"

[interior_energetics]
    module = "aragog"               # required: SPIDER has no re-melt path
```

Eight embryos between 0.05 and 0.11 AU at this spacing go unstable within a few thousand years, so the whole impact history is over long before the interior has cooled. That is deliberate: it makes the coupling visible in a short run. It is not a model of Earth's accretion, which plays out over 10<sup>7</sup> to 10<sup>8</sup> yr at 1 AU. Note also that these embryos need not orbit where the PROTEUS planet orbits; only the *fractional* orbit change is transferred, so the dynamical environment and the simulated planet's distance are independent choices.

During start-up, the run reports the system it evolved and the schedule it kept:

```
[ INFO  ] Running giant-impact model for 8 embryos
[ INFO  ] Following body 6 (selector 'mass'): 1.300 -> 2.500 M_earth, 0.0973 -> 0.0959 AU
[ INFO  ] Body 6 experienced 2 impacts
[ INFO  ] Scheduled 2 impact(s)
[ INFO  ]     first at 9.1183e+02 yr, last at 1.9394e+03 yr
```

Three of the eight embryos survive; `selector = "mass"` follows the heaviest, which grows from 1.3 to 2.5 M⊕ across two impacts, at about 912 yr (struck by body 5, adding 0.700 M⊕) and 1939 yr (struck by body 7, adding 0.500 M⊕).

As each impact lands, the loop reports it in a block that opens with the time-stepper announcing the shortened step and the impact itself (its target, impactor and added mass), and closes with the planet's resulting mass, orbital distance and eccentricity. Between those, one indented line is written per consequence that actually applied: the volatile and loss modes in force, the erosion fraction the law returned, what was stripped and what was delivered, the re-melt heat booked into the energy budget, and the mantle reset. How many of those appear depends on the configuration, so expect more lines under `match_planet` and `zephyrus` than under `dry` and `none`.

The added mass, the time and the two body ids come straight from the schedule above. The erosion fraction and the resulting planet state depend on the atmosphere the planet happens to have at that moment, so they differ from run to run even for the same schedule.

## Common pitfalls

**Aragog is required.** `interior_energetics.module = "spider"` with accretion enabled is refused at config load. SPIDER keeps its state in a restart file written by an external binary and has no validated re-melt path, so an impact cannot reset its mantle.

**No impacts scheduled.** The most common cause is `evolution_time` being short relative to the instability time of the configuration you chose. [Choosing the initial conditions](configuration.md#choosing-the-initial-conditions) covers how spacing and eccentricity set that time.

**Impacts before the run starts are discarded, not absorbed.** Morrigan measures time from disk dispersal. `time_offset` is added to each impact time, and any impact still landing at or before the run's start time is dropped with a warning naming the mass it would have added. That mass is **not** folded into the initial condition: the configured `planet.mass_tot` stands, and the planet ends the run lighter than the dynamics described. If the warning appears, shift `time_offset` so the history falls inside the simulated interval rather than accepting the truncated schedule.

**Spacing too wide for the masses.** Beyond a limit that depends on the embryo masses the layout is refused with an error naming the pair; see [choosing the initial conditions](configuration.md#choosing-the-initial-conditions).

**A settings file carrying `atm_mass_fraction`.** The model tracks no atmosphere; the key is ignored and a warning says so. Atmospheric loss belongs to `atmloss_module` on the PROTEUS side.

## Next step

- [Coupling to PROTEUS (explanation)](../Explanations/proteus_coupling.md) for what happens inside each of these steps.
- [Settings and inputs](../Reference/parameters.md) for the standalone settings file.
