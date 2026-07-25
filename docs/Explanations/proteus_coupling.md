# Coupling to PROTEUS

This page is the **theory** of how Morrigan plugs into a PROTEUS coupled run: where it sits in the loop, what the wrapper does with its output, how one planet is chosen out of a whole system, and what each impact changes when it lands. For the practical TOML recipe, see the [how-to page](../How-to/proteus_coupling.md).

## Morrigan runs once, not every iteration

Every other PROTEUS module is called inside the time loop and updates its slice of the planet's state on each pass. Morrigan is the exception. It runs **once**, before the loop starts, and what it returns is a **schedule**: a list of impacts with the times they occur and the physical parameters of each. The loop then consumes that schedule as it advances.

```
   ┌──────────────────────────────────────────────────────────────────────────┐
   │  init_accretion(handler)                    <-- Morrigan runs HERE, once │
   │      build_parameters(config)                   config + star.mass       │
   │      morrigan.run_system(**params)              evolve the whole system  │
   │      select_planet(survivors, config)           pick ONE survivor        │
   │      -> [ImpactEvent, ImpactEvent, ...]         its impact history       │
   │      validate_timeline(events)                  reject an impossible one │
   │      _drop_events_before_start(events, t0)      discard the pre-run tail │
   │                                                                          │
   │  while not done:                                                         │
   │                                                                          │
   │      pending = next_event(events, t)          # the loop reads the       │
   │      interior_o.t_next_impact = pending.time  # schedule, never the model│
   │                                                                          │
   │      run_interior(...)                        # cooling step; the        │
   │                                               # time-stepper shortens dt │
   │                                               # to land on t_next_impact │
   │                                                                          │
   │      hf_row['Time'] += dt                     # advance FIRST            │
   │                                                                          │
   │      for event in due_events(events, t_prev, t_now):   # impact falls due│
   │          apply_impact(handler, event)         # <-- consequences here    │
   │                                                                          │
   │      run_orbit(...) / run_escape(...) / run_outgassing(...) / ...        │
   └──────────────────────────────────────────────────────────────────────────┘
```

The time advance comes **before** the impact block, so the orbit, structure and
escape steps of the step an impact lands on already see the grown planet. The
clamp that shortens `dt` onto the impact time lives in the interior
time-stepper, not in the main loop; the loop only publishes when the next
impact is due.

The reason for the split is that the dynamical model and the coupled framework work on different clocks and at different cost. Morrigan integrates a whole system of embryos over hundreds of millions of years in a few seconds; PROTEUS integrates one planet's interior and atmosphere over the same span in hours. Running the dynamics once and replaying its schedule keeps the expensive loop in charge of the timestep.

A consequence worth stating plainly: **the coupling is one-way**. Nothing the coupled planet does feeds back into the dynamics. If the planet loses half its atmosphere to escape, the impact schedule does not change, because the schedule was fixed before the loop began.

## Selecting one planet out of a system

Morrigan evolves a system of embryos and typically leaves several survivors. PROTEUS models one planet. Something has to choose, and that choice is the coupling's main free decision.

`select_planet` implements four rules, set by `accretion.morrigan.selector`:

| Selector | Picks | Use when |
|---|---|---|
| `match_config` | the survivor whose **initial** mass and orbit are closest to the configured planet | you have a planet in mind and want the dynamics to describe it |
| `mass` | the survivor with the largest **final** mass | you want the system's dominant body |
| `semimajoraxis` | the survivor whose **final** orbit is nearest `selector_value`, in AU | you care about a particular orbital distance |
| `id` | the embryo with index `selector_value` | you are reproducing a specific system by hand |

Which end of the history a rule compares is what decides which body you get. `match_config` matches on the state the embryo *started* in, so it answers "which of these bodies began as the planet I configured"; `mass` and `semimajoraxis` match on the state it *ended* in, so they answer "which of these bodies became what I am looking for". A body can easily win one and lose the other.

The selected survivor's impact history becomes the schedule. Every other survivor is discarded. This is why `num_planets` and the initial mass distribution matter even though only one planet is followed: they set the dynamical environment that produced the one you keep.

## From Morrigan's record to a PROTEUS event

`build_parameters` translates the TOML block into the arguments `run_system` expects:

| `run_system` argument | Source | Note |
|---|---|---|
| `seed` | `accretion.morrigan.seed` | fixes the Monte Carlo |
| `masses` | `masses` or `num_planets` × `mass_equal` | converted M⊕ → kg |
| `eccentricity` | `eccentricity_init` | shared by all embryos |
| `inner_edge` | `inner_edge` | converted AU → m |
| `spacing` | `spacing` | mutual Hill radii |
| `density` | `density` | sets the mass-to-radius relation |
| `impact_angle` | `impact_angle` | its sine is the impact parameter |
| `evolution_time` | `evolution_time` | Gyr |
| `inner_cutoff` | `inner_cutoff` | converted AU → m |
| `stellar_mass` | **`star.mass`** | taken from the star block, not repeated here |

The stellar mass is deliberately not a Morrigan setting. It is read from `config.star.mass` so the dynamical model and the rest of PROTEUS cannot disagree about the host star.

Each returned record becomes an `ImpactEvent`. The one transformation applied on the way is the **time axis**: Morrigan measures time from disk dispersal, PROTEUS from the start of its own evolution, and `accretion.time_offset` maps between them. `time_offset` is added to each impact time, and the result is compared against the run's start time.

Impacts that still land at or before that start time are **discarded**, with a warning naming the mass they would have added. Their mass is not folded into the initial condition and arrives nowhere: the configured `planet.mass_tot` and orbit define the starting state on their own. A schedule whose early impacts fall outside the simulated interval therefore grows the planet less than the dynamics described, so choose `time_offset` to bring the history you care about inside the run.

`validate_timeline` then rejects a schedule that cannot describe one body: times must increase strictly, each impact's target mass must follow from the previous merged mass, and a body may not gain mass between impacts. A drop of up to 10 % between impacts is allowed, since a consumer may strip an atmosphere in between.

## What an impact does

When the loop reaches a scheduled time, `apply_impact` applies seven consequences in a fixed order. The order matters: each step reads state the previous one wrote.

1. **Size the volatile consequences** from the pre-impact state, before anything changes. The erosion fraction, what the target loses, and what the impactor carries are all computed first, so no step sees a half-updated planet.
2. **Grow the planet.** `planet.mass_tot` increases by the impactor's *rock* alone, the impactor's mass less its volatile content. The volatiles arrive through their own channel below, so counting them here would double-count them.
3. **Re-solve the structure.** The interior module recomputes radius, core size, and pressures at the new mass.
4. **Apply the volatile changes.** The target loses its stripped fraction, the impactor delivers what survives, and the whole-planet element budget is refreshed.
5. **Re-melt the mantle.** The interior is reset to its molten initial condition, recomputed for the grown planet, and the next interior solve is told not to clip the resulting temperature jump. On Aragog the heat this injects is measured across the cooled-to-molten entropy jump and booked into the energy budget as `step_dE_impact_J`, so an impact appears in the budget as a named source rather than as an unexplained residual. The impact's *kinetic* energy is logged alongside it for comparison but is not itself booked; the two are different quantities. How molten the reset state is depends on `planet.temperature_mode`: only `liquidus_super` solves for a profile that is guaranteed fully molten. `adiabatic_from_cmb` can be made molten by pinning `tcmb_init` high enough, and `accretion` usually is, but neither is checked against the liquidus. The remaining modes are only as molten as the configured initial condition, and the run warns about those at start-up.
6. **Clear the solidification latch.** A mantle that had crystallised is molten again, so the one-way latch is lifted. Without this step outgassing would stay frozen and the volatiles would be treated as locked in a solid mantle for good.
7. **Move the orbit.** The semi-major axis is scaled by the impact's *fractional* change and the eccentricity is set to its post-impact value.

The orbit step is worth expanding. Morrigan's absolute orbits belong to its own system, which need not sit where the PROTEUS planet sits. So the coupling applies the **ratio** `a_after / a_before`, not the absolute value: the planet keeps its configured orbital distance and inherits the dynamical model's fractional kick.

## What PROTEUS adds that Morrigan does not do

Morrigan reports bare bodies. Three pieces of physics live entirely on the PROTEUS side:

- **Impact atmospheric erosion.** Morrigan tracks no atmosphere and its merged masses are exact sums. PROTEUS computes what an impact strips, either as a fixed fraction or through the Kegerreis et al. (2020)[^cite-kegerreis2020] scaling law evaluated by `zephyrus.collision.mass_loss` from the record's own collision parameters. One fraction governs both bodies: the target loses that fraction of its atmosphere, and a volatile-bearing impactor loses the same fraction of its atmospheric part.
- **Volatile delivery.** Whether an impactor carries volatiles at all is a PROTEUS choice: dry, matching the planet's formation composition, or per-element budgets.
- **The thermal consequence.** Re-melting the mantle, booking the heat that injects, and lifting the solidification latch is interior physics, not dynamics.

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

`run_system` seeds numpy's global random state for its own draws and restores the caller's state afterwards, so an embedding program's random sequence is unaffected by running a system in-process. This matters in a coupled run, where PROTEUS and other modules draw from the same global state.

## See also

- [Coupling to PROTEUS (how-to)](../How-to/proteus_coupling.md) for the TOML recipe, the selector choice, and the pitfalls.
- [Model overview](model.md) for the dynamics that produce the schedule.
- [Limitations](limitations.md) for what the model does not represent.
- The PROTEUS-side code: `src/proteus/accretion/wrapper.py` (dispatch and impact consequences), `src/proteus/accretion/morrigan.py` (this module's adapter), `src/proteus/accretion/common.py` (the `ImpactEvent` schema and timeline validation), `src/proteus/config/_accretion.py` (the config block).

## References

[^cite-kegerreis2020]: Kegerreis, J.A., Eke, V.R., Catling, D.C., Massey, R.J., Teodoro, L.F.A. & Zahnle, K.J., *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020.
