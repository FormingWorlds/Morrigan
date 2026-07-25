# Installation

Morrigan is a pure-Python package with no compiled dependencies. Python 3.11 or newer on Linux or macOS is required; the `pylaplace` dependency, which supplies the Laplace coefficients for the secular solution, installs from PyPI.

## User install

```bash
git clone git@github.com:FormingWorlds/Morrigan.git
cd Morrigan
pip install -e .
```

## Developer install

The `develop` extra adds the test stack (pytest, coverage, hypothesis, the ZEPHYRUS cross-check dependency) and pre-commit:

```bash
pip install -e ".[develop]"
pre-commit install -f
```

## Optional extras

- `plot`: matplotlib and scipy for the repository's `plot.py` script, kept out of the base install so an embedding framework does not pull in a plotting stack.
- `docs`: the Zensical documentation toolchain, for building this site locally with `zensical serve`.

## Verifying the install

```bash
pytest -m "(unit or smoke) and not skip"
```

runs the fast test tiers, which need only the base install plus the `develop` extra; see [Testing](testing.md) for the full tier layout.
