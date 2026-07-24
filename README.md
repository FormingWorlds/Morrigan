# Morrigan

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/docs.yaml?branch=main&label=Docs)](https://github.com/FormingWorlds/Morrigan/actions/workflows/docs.yaml)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/Morrigan?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/Morrigan)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/Morrigan/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/Morrigan/actions/workflows/nightly.yml)

**Morrigan** is the dynamical evolution module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. This module contains the (modified) Python version of the semi-analytical model for giant impacts developed by [Kimura et al 2025](https://iopscience.iop.org/article/10.3847/1538-4357/ade992/meta) that models the orbital and accretionary evolution of planets as a result of giant impacts and gravitational scattering.

Named after Morrigan, a shapeshifting figure from Irish mythology thought to represent the dynamical nature of existence.

New feature(s) introduced in Python version:
1. Fractional atmospheric loss due to giant impact mergers

## Installation

```bash
git clone git@github.com:FormingWorlds/Morrigan.git
cd Morrigan
pip install -e .
```

## Running a model

Settings live in a `.toml` file; `initialise.toml` in the repository root is a worked example covering the number of systems, the embryo masses and atmospheric mass fractions, the stellar mass, and how long to evolve for.

```bash
morrigan -c initialise.toml
```

This replaces `python run_model.py`. If you have a launcher script of your own that still calls the old path, that is the line to change.

Keeping several settings files side by side and pointing at whichever you want is the reason for the `-c` flag. The repository ignores `*.toml` apart from `initialise.toml`, so a settings file you want to keep alongside the results it produced has to be added to git explicitly.

Results are written under the `save_directory` named in the settings file. A relative `save_directory` is taken from the directory you run the command in, not from wherever the settings file happens to live, so the run prints the full path it is writing to:

| Path | Contents |
| --- | --- |
| `data/full_systems/` | State of every planet through time |
| `data/mergers/` | One row per collision: the bodies involved, collision velocity, and atmospheric loss |
| `data/survivors/` | Final mass, orbit, eccentricity, and atmosphere of each surviving planet |
| `batch_summary.csv` | Runtime and surviving-planet count for each system |

To plot the results of a run, point the plotting script at the same settings file:

```bash
pip install -e ".[plot]"
python plot.py -c initialise.toml
```

`plot.py` is a script kept in the repository rather than part of the installed package, so it is run from a checkout. Its plotting dependencies (matplotlib, scipy) are in the `plot` extra, which keeps them out of the way of anything that only wants to import the model.

## Documentation

The full documentation, including the model physics, the PROTEUS coupling contract, and the per-source validation anchors, lives at [proteus-framework.org/Morrigan](https://proteus-framework.org/Morrigan). Build it locally with `pip install -e ".[docs]"` and `zensical serve`.

## Testing

```bash
pip install -e ".[develop]"
pytest -m "(unit or smoke) and not skip"
```

The suite is tiered (unit and smoke on every pull request, seed-ensemble statistics nightly) and every physics routine is pinned against a published value, an analytical limit, or a cross-implementation check; see the documentation's testing guide for the full contract.

## Reproducibility

Each system is seeded from `random_seed` in the settings file plus its own index, so a given settings file reproduces the same systems exactly.
