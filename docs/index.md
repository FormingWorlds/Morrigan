# Morrigan

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/docs.yaml?branch=main&label=Docs)](https://github.com/FormingWorlds/Morrigan/actions/workflows/docs.yaml)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/Morrigan?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/Morrigan)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/nightly.yml)

**Morrigan** is the giant-impact accretion module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It implements the semi-analytical Monte Carlo model of [Kimura et al. (2025)](https://doi.org/10.3847/1538-4357/ade992) for the post-disk dynamical evolution of planetary systems: a set of planetary embryos evolves through secular eccentricity oscillations, orbit crossings, gravitational scatterings, ejections, and giant impacts, until the system settles into a Hill-stable configuration.

Morrigan runs standalone from a TOML settings file, producing per-system tables of the full evolution, every merger, and the surviving planets, and it runs in memory through `morrigan.run_system`, the entry point PROTEUS uses to drive its accretion coupling: each reported impact grows the coupled planet, re-melts its mantle, delivers volatiles, erodes its atmosphere, and moves its orbit.

Named after Morrigan, a shapeshifting figure from Irish mythology thought to represent the dynamical nature of existence.

## Where to go

- [Installation](How-to/installation.md) and [running a model](How-to/running.md)
- [Model overview](Explanations/model.md): the physics, equation by equation
- [Coupling to PROTEUS](Explanations/proteus_coupling.md): the impact-record interface
- [Settings and inputs](Reference/parameters.md) and the [API reference](Reference/api.md)
- [Validation anchors](Validation/index.md): where every physics routine is pinned against a published value, an analytical limit, or a cross-implementation check

## The model in one paragraph

Embryos are laid out with a fixed spacing in mutual Hill radii and evolve under the classical Laplace-Lagrange secular solution. The time to the next orbit crossing is predicted from the Petit et al. (2020) three-body resonance-overlap timescale, gated by a two-planet stability criterion. When a crossing falls due, the interacting pair either collides, with a probability built from the ratio of its relative eccentricity to the mutual escape eccentricity, or scatters; scattering can excite an orbit past unity eccentricity and eject a body. Collisions merge the pair perfectly, summing the masses and placing the merged body on the orbital-energy-weighted orbit; the model tracks no atmosphere, so impact erosion is left to the consumer. The run ends when the system is Hill-stable or the evolution time is spent.
