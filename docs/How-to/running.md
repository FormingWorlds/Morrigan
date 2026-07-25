# Running a model

## From the command line

Settings live in a TOML file; `initialise.toml` in the repository root is a worked example. All settings are documented in [Settings and inputs](../Reference/parameters.md).

```bash
morrigan -c initialise.toml
```

Keeping several settings files side by side and pointing at whichever you want is the reason for the `-c` flag. The repository ignores `*.toml` apart from `initialise.toml`, so a settings file you want to keep alongside the results it produced has to be added to git explicitly.

Results are written under the `save_directory` named in the settings file. A relative `save_directory` is taken from the directory you run the command in, not from wherever the settings file lives, so the run prints the full path it is writing to:

| Path | Contents |
| --- | --- |
| `data/full_systems/` | State of every planet through time |
| `data/mergers/` | One row per collision: the bodies involved and the collision velocity |
| `data/survivors/` | Final mass, orbit, and eccentricity of each surviving planet |
| `batch_summary.csv` | Runtime and surviving-planet count for each system |

A batch of several systems (the `[batch]` table in the settings file) runs on a process pool, one system per worker, each seeded from the base seed plus its own index.

## Plotting

```bash
pip install -e ".[plot]"
python plot.py -c initialise.toml
```

`plot.py` is a script kept in the repository rather than part of the installed package, so it is run from a checkout.

## In memory, from another program

`morrigan.run_system` evolves one system without writing files and returns the survivors and per-impact records; this is the interface PROTEUS drives. See [Coupling to PROTEUS](../Explanations/proteus_coupling.md) for the call signature and the record schema.

## Reproducibility

Each system is seeded from `random_seed` in the settings file plus its own index, so a given settings file reproduces the same systems exactly. The same holds for `run_system` with a fixed seed: one seed, one timeline, bit for bit.
