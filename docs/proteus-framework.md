---
title: Morrigan in the PROTEUS framework
---

<h1 align="center">
    <a href="https://proteus-framework.org">
    <div>
        <img src="https://raw.githubusercontent.com/FormingWorlds/PROTEUS/main/docs/assets/PROTEUS_white.png#gh-light-mode-only" style="vertical-align: middle;" width="60%"/>
        <img src="https://raw.githubusercontent.com/FormingWorlds/PROTEUS/main/docs/assets/PROTEUS_black_nobkg.png#gh-dark-mode-only" style="vertical-align: middle;" width="60%"/>
    </div>
    </a>
</h1>


Morrigan is the **protoplanet accretion module** of [PROTEUS](https://proteus-framework.org/PROTEUS) (/ˈproʊtiəs/, PROH-tee-əs), a modular Python framework for the coupled evolution of the atmospheres and interiors of rocky planets and exoplanets. A schematic of PROTEUS components and corresponding modules can be found below. Click any module in the diagram to open its documentation, or navigate to it from the sidebar.
<br>

<object type="image/svg+xml" data="https://cdn.jsdelivr.net/gh/FormingWorlds/PROTEUS@main/docs/assets/proteus_modules_schematic.svg" class="mod-diagram mod-diagram--light" aria-label="PROTEUS module schematic (light mode)">PROTEUS module schematic (light mode)</object>
<object type="image/svg+xml" data="https://cdn.jsdelivr.net/gh/FormingWorlds/PROTEUS@main/docs/assets/proteus_modules_schematic_darkmode.svg" class="mod-diagram mod-diagram--dark" aria-label="PROTEUS module schematic (dark mode)">PROTEUS module schematic (dark mode)</object>

<p style="text-align: center;"><strong>Schematic of PROTEUS components and corresponding modules.</strong></p>

## Where Morrigan sits


Morrigan occupies the **accretion** slot. It answers the question the rest of the loop cannot: what did this planet's assembly look like, and when did it get hit?

Unlike the other modules it is not called on each iteration. It runs once, before the loop starts, and returns a schedule of impacts that the loop then replays. See [Coupling to PROTEUS](Explanations/proteus_coupling.md) for why.

## What Morrigan exchanges with

The schematic above links every module. Morrigan interacts directly with two of them:

| Module | Exchange |
| --- | --- |
| [MORS](https://proteus-framework.org/MORS/) and the stellar block | Morrigan reads the host star's mass from `star.mass` rather than repeating it, so the dynamics and the rest of the run cannot disagree about the star |
| [ZEPHYRUS](https://proteus-framework.org/ZEPHYRUS/) | Morrigan models no atmosphere, so impact erosion is left to the consumer. When a coupled run selects `atmloss_module = "zephyrus"`, PROTEUS passes each scheduled impact's collision parameters to ZEPHYRUS's collision law to work out how much atmosphere it strips. The call is PROTEUS's; Morrigan has no ZEPHYRUS dependency |

## Conventions the ecosystem shares

- **SI internally.** Modules exchange quantities in kg, m and s. Configuration files use the units a user thinks in, converted once at the boundary.
- **One module per slot per run.** Where several modules serve the same slot they are alternatives, chosen in the configuration.
- **The interface is the contract.** A module's public entry point and the schema it returns are what the framework depends on; changing either is a breaking change for the coupled run and is flagged as such.
