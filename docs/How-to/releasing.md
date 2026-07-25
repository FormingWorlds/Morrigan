# Releasing

Morrigan's version comes from git tags through setuptools-scm, so tagging is what makes a release. There is no version string to edit anywhere in the source.

## Version scheme

CalVer, `YY.0M.0D`, matching the rest of the PROTEUS ecosystem. A release cut on 3 August 2026 is tagged `26.08.03`.

The scheme carries no compatibility promise in the numbers themselves. What changes in a way that breaks a consumer is the impact-record schema, and that is called out in the pull request that changes it rather than encoded in the version. PROTEUS pins the Morrigan it builds against, so a release that changes the schema needs that pin moved in the same change.

## Before tagging

1. `main` is green: the test workflow passes on the commit you intend to tag.
2. The documentation builds without warnings, since the docs deploy runs on push to `main`.
3. Anything that changed the record schema, the units, or the physics is described in its pull request, because the release notes are assembled from them.

## Cutting the release

```bash
git checkout main
git pull
git tag 26.08.03
git push origin 26.08.03
```

Then create the GitHub release against that tag, either in the web interface or with

```bash
gh release create 26.08.03 --repo FormingWorlds/Morrigan \
    --title 26.08.03 --generate-notes
```

`--generate-notes` assembles the notes from the merged pull requests, which is why their descriptions matter.

## Verifying

```bash
pip install -e .
python -c "import morrigan; print(morrigan.__version__)"
```

should print the tag you just pushed. A checkout with no tags in its history reports `0.0.0`, the configured fallback, which is the usual sign that a shallow clone dropped the tag history rather than that anything is wrong with the release.

## Publication

Publishing is automatic. Creating the GitHub release triggers the publish workflow, which builds the distribution and uploads it to PyPI as `fwl-morrigan` through trusted publishing, so no API token is stored in the repository.

The workflow refuses to upload a fallback version. setuptools-scm reports `0.0.0` when it cannot see the tag history, and a shallow checkout is the usual cause; publishing that would take the name `0.0.0` on PyPI permanently, and PyPI does not allow a version to be re-uploaded. The build is checked before the upload step for that reason.

Confirm the upload landed by querying the index rather than by installing:

```bash
curl -s https://pypi.org/pypi/fwl-morrigan/json | python -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
```

Installing to check would replace an editable development checkout with the released wheel, which is rarely what you want mid-session.

Note that PyPI normalises the CalVer tag: `26.08.03` is published as `26.8.3`. Both refer to the same release, and a version-based dependency pin should use the normalised form.

## Moving PROTEUS onto a new release

PROTEUS records the Morrigan it builds against in its own `pyproject.toml`, under `[tool.proteus.modules.morrigan]`, as a repository URL plus a `ref` naming an exact commit. `tools/get_morrigan.sh` resolves that entry and installs the checkout editable, so the `ref` is what decides which Morrigan a coupled run uses. After a release, update it there and open a pull request against PROTEUS.

A commit `ref` pins a specific tree rather than a version range, so it does not move on its own; a release that is not followed by that update leaves coupled runs on the previous commit.
