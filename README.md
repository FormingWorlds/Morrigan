# Morrigan

**Morrigan** is the dynamical evolution module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. This module contains the (modified) Python version of the semi-analytical model for giant impacts developed by [Kimuara et al 2025](https://iopscience.iop.org/article/10.3847/1538-4357/ade992/meta) that models the orbital and accretionary evolution of planets as a result of giant impacts and gravitational scattering.

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

Results are written under the `save_directory` named in the settings file:

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

## Reproducibility

Each system is seeded from `random_seed` in the settings file plus its own index, so a given settings file reproduces the same systems exactly.
