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

## Worked example: an accreting Earth analogue

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

A run like this reports the schedule it built during start-up, then each impact as it lands:

```
[ INFO  ] Running giant-impact model for 8 embryos
[ INFO  ] Body 6 experienced 1 impacts
[ INFO  ] Scheduled 1 impact(s)
...
[ INFO  ] Time-stepping: impact at 9.5991e+02 yr, capping dt at 1.00e+02 yr
[ INFO  ] Giant impact at t = 9.5991e+02 yr: target 6 struck by 7, adding 0.5000 M_earth
[ INFO  ]     impactor volatiles: match_planet; atmosphere loss: zephyrus
[ INFO  ]     impact erosion law: loss fraction 0.554
[ INFO  ]     planet is now 1.4999 M_earth at 1.03107 AU, e = 0.0223
```

## Common pitfalls

**Aragog is required.** `interior_energetics.module = "spider"` with accretion enabled is refused at config load. SPIDER keeps its state in a restart file written by an external binary and has no validated re-melt path, so an impact cannot reset its mantle.

**No impacts scheduled.** The most common cause is `evolution_time` being short relative to the instability time of the configuration you chose. [Choosing the initial conditions](configuration.md#choosing-the-initial-conditions) covers how spacing and eccentricity set that time.

**All impacts before the run starts.** Morrigan measures time from disk dispersal. If every impact lands before your `time_offset`, they are folded into the initial condition and the run proceeds with none scheduled. Check the reported impact times against your start time.

**Spacing too wide for the masses.** Beyond a limit that depends on the embryo masses the layout is refused with an error naming the pair; see [choosing the initial conditions](configuration.md#choosing-the-initial-conditions).

**A settings file carrying `atm_mass_fraction`.** The model tracks no atmosphere; the key is ignored and a warning says so. Atmospheric loss belongs to `atmloss_module` on the PROTEUS side.

## Next step

- [Coupling to PROTEUS (explanation)](../Explanations/proteus_coupling.md) for what happens inside each of these steps.
- [Settings and inputs](../Reference/parameters.md) for the standalone settings file.
