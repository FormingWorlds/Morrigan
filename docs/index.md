# Morrigan

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/docs.yaml?branch=main&label=Docs)](https://github.com/FormingWorlds/Morrigan/actions/workflows/docs.yaml)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/Morrigan?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/Morrigan)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/nightly.yml)

**Morrigan** is the protoplanet accretion module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It follows a system of protoplanets through the giant impacts and gravitational scattering by which they accrete, implementing the semi-analytical Monte Carlo model of [Kimura et al. (2025)](https://doi.org/10.3847/1538-4357/ade992). The protoplanets are the planetary embryos of the post-disk stage, and that is what the model and this documentation call them: a set of embryos evolves through secular eccentricity oscillations, orbit crossings, gravitational scatterings, ejections, and giant impacts, until the system settles into a Hill-stable configuration.

Morrigan runs standalone from a TOML settings file, producing per-system tables of the full evolution, every merger, and the surviving planets, and it runs in memory through `morrigan.run_system`, the entry point PROTEUS uses to drive its accretion coupling: each reported impact grows the coupled planet, re-melts its mantle, delivers volatiles, erodes its atmosphere, and moves its orbit.

Named after Morrigan, a shapeshifting figure from Irish mythology thought to represent the dynamical nature of existence.

## Where to go

<div class="grid cards" markdown>

-   :material-rocket-launch: **Get started**

    The shortest path from a clone to a result

    [Go to getting started](getting_started.md)

-   :material-school: **Follow a tutorial**

    One system standalone, or one inside PROTEUS

    [Go to the tutorials](Tutorials/standalone.md)

-   :material-download: **Install**

    From PyPI, or as an editable checkout

    [Go to installation guide](How-to/installation.md)

-   :material-play: **Run a model**

    Drive a system from a settings file

    [Go to running a model](How-to/running.md)

-   :material-tune: **Configure**

    Choose the initial conditions and what they do

    [Go to configuration](How-to/configuration.md)

-   :material-book-open-variant: **Understand the model**

    The physics, equation by equation

    [Go to model overview](Explanations/model.md)

-   :material-sitemap: **Find it in the source**

    Where each piece of the model lives

    [Go to code architecture](Explanations/code_architecture.md)

-   :material-link-variant: **Couple to PROTEUS**

    The TOML recipe and the pitfalls

    [Go to the how-to](How-to/proteus_coupling.md)

-   :material-lightbulb-on: **Understand the coupling**

    Why it is built the way it is

    [Go to the explanation](Explanations/proteus_coupling.md)

-   :material-table: **Look up a setting**

    Every field of the settings file

    [Go to settings and inputs](Reference/parameters.md)

-   :material-code-braces: **Browse the API**

    The entry points a caller uses

    [Go to API reference](Reference/api.md)

-   :material-check-decagram: **See what is validated**

    Each routine pinned against a published or analytical value

    [Go to validation anchors](Validation/index.md)

-   :material-alert-circle-outline: **Know the limits**

    What the model does not represent

    [Go to limitations](Explanations/limitations.md)

-   :material-github: **Browse the code**

    The repository and its history

    [Go to source code](https://github.com/FormingWorlds/Morrigan)

-   :material-bug: **Raise an issue**

    Report a problem or ask for a feature

    [Go to issues](https://github.com/FormingWorlds/Morrigan/issues)

</div>

## The model in one paragraph

Embryos are laid out with a fixed spacing in mutual Hill radii and evolve under the classical Laplace-Lagrange secular solution. The time to the next orbit crossing is predicted from the Petit et al. (2020) three-body resonance-overlap timescale, gated by a two-planet stability criterion. When a crossing falls due, the interacting pair either collides, with a probability built from the ratio of its relative eccentricity to the mutual escape eccentricity, or scatters; scattering can excite an orbit past unity eccentricity and eject a body. Collisions merge the pair perfectly, summing the masses and placing the merged body on the orbital-energy-weighted orbit; the model tracks no atmosphere, so impact erosion is left to the consumer. The run ends when the system is Hill-stable or the evolution time is spent.
