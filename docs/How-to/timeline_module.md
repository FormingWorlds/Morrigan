# Replay an impact history: the `timeline` module

`accretion.module = "timeline"` applies a sequence of impacts read from a file
instead of deriving one from a dynamical model. Every consequence is computed
exactly as it is for a Morrigan-derived history, so the two paths differ only in
where the schedule came from.

This is a PROTEUS module, not part of Morrigan. It is documented here because a
Morrigan run produces exactly the file it consumes, and because it is how a
Morrigan-derived history reaches someone who does not have Morrigan installed.

## When to reach for it

- **Reproduce a published history**, or one a collaborator sends you, without
  rerunning the dynamics.
- **Share a run** with someone who does not have the optional Morrigan
  dependency: the file is all they need.
- **Isolate one impact.** A hand-written two-line file is the cleanest way to
  ask what a single collision does to an interior, with no Monte Carlo in the
  way.
- **Hold the schedule fixed** while varying something else, so two runs differ
  only in the quantity under study.

## Configuration

```toml
[accretion]
    module      = "timeline"
    time_offset = 0.0

    [accretion.timeline]
        timeline_path = "output/my_run/impact_timeline.csv"
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `timeline_path` | str or none | `none` | Path to the impact timeline file. Environment variables and `~` are expanded |

`timeline_path` is required for this module and refused on the others, so a
configuration that names a file but forgets to select the module is rejected at
load rather than quietly running something else.

## Where the file comes from

Every PROTEUS run with accretion enabled writes the schedule it resolved to
`impact_timeline.csv` in its output directory, whichever module produced it.
That file is already on the PROTEUS time axis, with `accretion.time_offset`
applied when it was written.

!!! warning "Do not offset it twice"
    `time_offset` is added again when a timeline file is read. Leave it at its
    default of `0.0` when replaying a file that PROTEUS wrote, or the whole
    history shifts by the offset a second time.

## The file format

Comma- or whitespace-separated, one header row, `#` starts a comment. Every
column below must be present:

| Column | Unit | Meaning |
|---|---|---|
| `time` | yr | when the impact occurs |
| `M_target_before` | kg | target mass before the collision |
| `M_impactor` | kg | impactor mass |
| `M_merged_after` | kg | perfect-merger sum of the two |
| `v_impact` | m s<sup>-1</sup> | speed at first contact |
| `v_esc` | m s<sup>-1</sup> | mutual escape speed |
| `impact_parameter` | 1 | sine of the impact angle; 0 head-on, 1 grazing |
| `R_target_before` | m | target radius before the collision |
| `R_impactor` | m | impactor radius |
| `rho_target` | kg m<sup>-3</sup> | target bulk density |
| `rho_impactor` | kg m<sup>-3</sup> | impactor bulk density |
| `a_before` | m | semi-major axis before |
| `a_after` | m | semi-major axis after |
| `e_before` | 1 | eccentricity before |
| `e_after` | 1 | eccentricity after |
| `id_target` | 1 | target body id |
| `id_impactor` | 1 | impactor body id |

The bodies carry no atmosphere: masses and radii are bulk values, and the
densities follow from them, so each row is self-consistent on its own.

Both orbital elements are recorded before and after, because the coupling
transfers the *change* rather than the absolute orbit: the ratio for the
semi-major axis, the difference for the eccentricity. A file that reports only
the post-impact values cannot drive the coupling.

## What is checked when it loads

A timeline has to describe one body's history, and a file that cannot is
rejected rather than run:

- times increase strictly;
- each impact's `M_target_before` follows from the previous row's
  `M_merged_after`;
- a body does not gain mass between impacts. A drop of up to 10 % is allowed,
  since a consumer may strip an atmosphere in between.

Impacts landing at or before the run's start time are discarded with a warning
naming the mass they would have added. That mass is not folded into the initial
condition: the configured `planet.mass_tot` stands.

## See also

- [Coupling to PROTEUS (how-to)](proteus_coupling.md) for the recipe a
  Morrigan-driven run uses.
- [The analytical `dummy` module](dummy_module.md) for a history from scaling
  laws rather than from a file.
- [Coupling to PROTEUS (explanation)](../Explanations/proteus_coupling.md) for
  what each impact does when it lands.
