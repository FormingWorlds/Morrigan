# Getting started

Morrigan can be used two ways: as a standalone model that writes tables of a system's evolution, or in-process as the accretion module of a PROTEUS run. Both come from the same install.

## Quick path

```bash
git clone git@github.com:FormingWorlds/Morrigan.git
cd Morrigan
pip install -e .
morrigan -c initialise.toml
```

That evolves the worked example in `initialise.toml` and writes its results under the `save_directory` named there. The run prints the full path it is writing to.

## What do you want to do?

**Run a system and look at the results.** Start with [running a model](How-to/running.md), which covers the settings file, the output tables and the plotting script. [Settings and inputs](Reference/parameters.md) is the full field reference.

**Understand what the model does.** [Model overview](Explanations/model.md) walks the physics in the order the code applies it: how embryos are laid out, how the next instability is timed, and how a crossing resolves into a collision, a scattering or an ejection. [Limitations](Explanations/limitations.md) is the companion: what the model deliberately does not represent.

**Drive a PROTEUS run with it.** [Coupling to PROTEUS (how-to)](How-to/proteus_coupling.md) is the TOML recipe. [Coupling to PROTEUS (explanation)](Explanations/proteus_coupling.md) is why it is built that way, including the one property that makes Morrigan unlike every other PROTEUS module: it runs once, before the loop, and produces a schedule rather than updating state each iteration.

**Change the model.** [Code architecture](Explanations/code_architecture.md) maps the source files onto the physics. [Testing](How-to/testing.md) is the contract any change has to meet, and [Validation anchors](Validation/index.md) records what each physics routine is pinned against. [Contributing](Community/CONTRIBUTING.md) covers what a change needs before it can land, including the two extra requirements on a physics change.

## Where the numbers come from

The model is the semi-analytical Monte Carlo of Kimura et al. (2025), with the ejection prescription from its companion application paper. Rather than integrating orbits, it predicts *when* a system becomes unstable and resolves each instability with analytic prescriptions calibrated against N-body simulations. That makes a single system cheap enough to run inside a coupled framework, at the cost of any individual history being a statistical realisation rather than a trajectory. Ensemble statistics are the meaningful comparison.
