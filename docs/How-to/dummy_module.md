# Grow a planet without dynamics: the `dummy` module

`accretion.module = "dummy"` builds an impact history from scaling laws instead
of integrating a system of embryos. It is deterministic: the same configuration
always produces the same history, with no Monte Carlo and no optional
dependency.

This is a PROTEUS module, not part of Morrigan. It is documented here because it
occupies the same slot, answers the same interface, and is usually the right
tool when what you want is a growing planet rather than a dynamical history.

## When to reach for it

- **Exercise the coupling** in a test or a quick run without installing
  Morrigan.
- **Prescribe the growth** you want to study, rather than sampling whatever a
  dynamical model produces.
- **Sweep an accretion parameter** where a single controlled history is worth
  more than the spread over seeds.

If what you need is a plausible *dynamical* history, with the scatterings,
ejections and the stochasticity that come with it, use Morrigan.

## Configuration

```toml
[accretion]
    module = "dummy"

    [accretion.dummy]
        mass_accreted    = 0.1       # [M_earth], total delivered over the timeline
        num_impacts      = 3
        timescale        = 1.0e6     # [yr], e-folding time of the growth law
        time_last        = 5.0e6     # [yr], time of the final impact
        eccentricity     = 0.05
        impact_parameter = 0.5
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `mass_accreted` | float | `0.1` | Total mass delivered over the whole timeline \[M<sub>&#8853;</sub>]; the increments are scaled to sum to exactly this |
| `num_impacts` | int | `3` | Number of impacts in the timeline |
| `timescale` | float | `1.0e6` | E-folding time of the accretion law \[yr] |
| `time_last` | float | `5.0e6` | Time of the final impact \[yr]. Impacts are spaced evenly from `time_last / num_impacts` up to this time |
| `eccentricity` | float | `0.05` | Encounter eccentricity \[1], setting the approach velocity added to the mutual escape velocity, and the impactor's orbit |
| `impact_parameter` | float | `0.5` | Impact parameter of every collision \[1], the sine of the impact angle. Zero is head-on, one is grazing |

## What the law does

The planet approaches an asymptotic mass exponentially, which is the standard
picture of an accretion rate that decays as the feeding zone empties. Impacts
are placed at evenly spaced times, and each delivers the mass the law accretes
over its interval, so **the increments decay with time and the first impact is
the largest**. The increments are then scaled so they sum to `mass_accreted`
exactly.

Everything else follows from the pair: radii from the Noack &amp; Lasbleis (2020)
mass-radius scaling, collision velocities from the mutual escape velocity
combined with an encounter velocity set by `eccentricity`, and the merged orbit
from conserving linear momentum.

`timescale` against `time_last` is what shapes the history. Short compared with
`time_last` concentrates almost all the mass in the first impact; long compared
with it spreads the mass nearly evenly across the impacts.

## Traps

!!! warning "`timeline_path` belongs to the other module"
    Setting `accretion.dummy.timeline_path` is refused at configuration load.
    The field exists only to catch that mistake: it asks to replay a file, and
    without the check the run would quietly generate a timeline at default
    settings instead. To replay a file, use
    [`accretion.module = "timeline"`](timeline_module.md).

The impacts are still real impacts as far as the rest of PROTEUS is concerned:
they grow the planet, deliver volatiles under `impactor_volatiles`, strip
atmosphere under `atmloss_module`, re-melt the mantle and move the orbit. Only
the schedule is analytical.

## See also

- [Replay an impact history](timeline_module.md) for a history read from a file.
- [Coupling to PROTEUS (how-to)](proteus_coupling.md) for the settings shared by
  all three modules.
- [Run PROTEUS with Morrigan](../Tutorials/proteus_minimal.md) for a complete
  configuration file.
