# Contributing guidelines

Morrigan is part of the [PROTEUS](https://proteus-framework.org/PROTEUS) ecosystem and follows its conventions. This page covers what a change has to satisfy before it can land.

## Getting set up

```console
git clone git@github.com:FormingWorlds/Morrigan.git
cd Morrigan
pip install -e ".[develop]"
pre-commit install -f
```

## Building the documentation

The documentation is written in [Markdown](https://www.markdownguide.org/basic-syntax/) and built with [Zensical](https://zensical.org/), which reads the existing `mkdocs.yml` directly.

```console
pip install -e ".[docs]"
zensical serve
```

For a one-shot production build:

```console
zensical build --clean
```

The build must finish with no warnings. A warning usually means a link that does not resolve or a footnote that is never referenced, and both render as visible defects on the published site.

## What a change has to satisfy

```console
pytest -m "(unit or smoke) and not skip"   # the tiers that run on every pull request
ruff check --fix src/ tests/ tools/
ruff format src/ tests/ tools/
bash tools/validate_test_structure.sh      # every test carries a tier marker
python tools/check_test_quality.py --check  # anti-happy-path lint
```

All four run in CI and block a merge. The test contract itself is documented in [Testing](https://proteus-framework.org/Morrigan/How-to/testing.html), and it is stricter than most: every physics routine must be pinned against a published value, an analytical limit or a cross-implementation check, and a pinned value needs a discrimination guard showing that the plausible wrong formula would fail the test.

A pull request that adds or substantially changes more than fifty lines of test code triggers an independent review before it can merge.

## Changing the physics

Two things beyond the tests:

**Cite the equation.** A change that transcribes or corrects a formula names the equation number in the source paper, in the pull request and in a comment at the call site. Transcription slips are this model's most common defect class, and the citation is what makes the next reader able to check it.

**Say what moved.** Most physics changes shift ensemble outcomes. State what you measured, over how many seeds, and in which configuration, so that someone comparing results across versions knows what to expect.

## Changing the record schema

The per-impact record is the interface PROTEUS consumes. Any change to its fields, units or conventions is a breaking change for a coupled run and has to be called out in the pull request that makes it, so the PROTEUS-side pin can be updated in step. See [Coupling to PROTEUS](https://proteus-framework.org/Morrigan/Explanations/proteus_coupling.html).

## Writing style

Documentation, comments and commit messages describe the current state of the code, not how it came to be. A reader arriving cold needs to know what is true now; the history is in the commit log where it belongs.

## Reporting a problem

Open an issue at [FormingWorlds/Morrigan/issues](https://github.com/FormingWorlds/Morrigan/issues). For a physics problem, the most useful report contains the settings file, the seed, and what you expected instead, because a single system is a statistical realisation and the seed is what makes it reproducible.
