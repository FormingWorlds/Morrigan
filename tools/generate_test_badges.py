"""Generate shields.io endpoint-badge JSON files for Morrigan test counts.

The script invokes ``pytest --collect-only -q`` per marker expression to
count tests without executing them, then writes one JSON file per badge
under the ``--out`` directory in the shields.io endpoint-badge schema:

    {"schemaVersion": 1, "label": "<text>", "message": "<count>", "color": "blue"}

Public-surface output (three files, two categories):

- ``tests-total.json`` (label "tests"): count of ``not skip`` tests.
- ``tests-unit.json`` (label "unit tests"): count of ``unit and not skip``.
- ``tests-integration.json`` (label "integration tests"): count of
  ``(smoke or integration or slow) and not skip``.

The pytest marker scheme registered in ``pyproject.toml`` has four tiers
(``unit``, ``smoke``, ``integration``, ``slow``); the public badge surface
intentionally collapses ``smoke + integration + slow`` into a single
"integration tests" category because a four-way taxonomy is confusing to
non-developer readers. The four markers remain available as per-test
decorators regardless of this collapse.

Every Morrigan test carries exactly one tier marker: ``unit``, ``smoke``,
``integration`` or ``slow``. ``tools/validate_test_structure.sh`` enforces
that and runs as a blocking check on every pull request, so the two
sub-categories partition the total and ``unit + integration == tests``
always holds. ``main`` asserts it, which turns the drift a newly
registered fifth tier would introduce into a loud failure instead of a
headline count that quietly stops matching the two beneath it.

Counting is by collection, so a test that pytest collects but skips at run
time (``importorskip``, ``skipif``) is still counted. Only the ``skip``
marker is filtered.

Usage
-----
    python tools/generate_test_badges.py --out badge_payload/

Notes
-----
Running the script does not execute the test suite; only collection is
triggered. Pytest exit code 5 ("no tests collected") is treated as a
successful zero count and the corresponding badge writes ``"message": "0"``.
Any other non-zero exit is a hard failure. Collection only reflects tests
whose modules import successfully, so the package must be installed (with
its ``develop`` extra) before this script runs.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

_COLLECT_RE = re.compile(r'^(\d+)(?:/\d+)?\s+tests?\s+collected\b', re.MULTILINE)

_BADGES: tuple[tuple[str, str, str], ...] = (
    ('total', 'tests', 'not skip'),
    ('unit', 'unit tests', 'unit and not skip'),
    (
        'integration',
        'integration tests',
        '(smoke or integration or slow) and not skip',
    ),
)

# Filenames not in the public-surface set; removed from the output
# directory at the end of every run so a reused output directory stays in
# sync with the three-file scheme above. The workflow writes to a fresh
# directory each run, so this only matters when running the tool by hand.
_PRUNE_FILES: tuple[str, ...] = ('tests-smoke.json', 'tests-slow.json')


def count_tests(marker_expr: str) -> int:
    """Run pytest collection and return the number of selected tests.

    Parameters
    ----------
    marker_expr : str
        Pytest marker expression passed via ``-m``.

    Returns
    -------
    int
        Number of tests pytest collected for the given marker. Exit
        code 5 ("no tests collected") is mapped to 0.

    Raises
    ------
    RuntimeError
        If pytest exits with a non-zero code other than 5, or if the
        trailing summary line cannot be parsed from stdout.
    """
    proc = subprocess.run(
        [sys.executable, '-m', 'pytest', '--collect-only', '-q', '-m', marker_expr],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode == 5:
        return 0
    if proc.returncode != 0:
        raise RuntimeError(
            f'pytest --collect-only -m {marker_expr!r} exited with '
            f'code {proc.returncode}\n'
            f'--- stdout ---\n{proc.stdout}\n'
            f'--- stderr ---\n{proc.stderr}'
        )
    match = _COLLECT_RE.search(proc.stdout)
    if match is None:
        raise RuntimeError(
            f'pytest summary line not found for marker {marker_expr!r}\n'
            f'--- stdout ---\n{proc.stdout}'
        )
    return int(match.group(1))


def write_badge(out_dir: Path, name: str, label: str, count: int) -> Path:
    """Write a shields.io endpoint-badge JSON file.

    Parameters
    ----------
    out_dir : Path
        Directory to write the JSON file into. Must already exist.
    name : str
        Suffix used in the filename ``tests-<name>.json``.
    label : str
        Badge label rendered on the left side of the shield.
    count : int
        Badge message count rendered on the right side of the shield.

    Returns
    -------
    Path
        Path of the written JSON file.
    """
    payload = {
        'schemaVersion': 1,
        'label': label,
        'message': str(count),
        'color': 'blue',
    }
    out_path = out_dir / f'tests-{name}.json'
    out_path.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    return out_path


def prune_extra_files(out_dir: Path) -> list[str]:
    """Remove any badge JSON files outside the public-surface set.

    Parameters
    ----------
    out_dir : Path
        Directory the JSON files live in.

    Returns
    -------
    list[str]
        Names of files that were removed, in the order encountered.
    """
    removed: list[str] = []
    for name in _PRUNE_FILES:
        path = out_dir / name
        if path.exists():
            path.unlink()
            removed.append(name)
    return removed


def main() -> int:
    """Entry point.

    Returns
    -------
    int
        Process exit code (always 0 on success; failures raise).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Directory to write tests-<name>.json badge files into.',
    )
    args = parser.parse_args()
    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    counts: dict[str, int] = {}
    for name, label, expr in _BADGES:
        count = count_tests(expr)
        counts[name] = count
        write_badge(out_dir, name, label, count)
        print(f'{label}: {count}')

    # Every test carries exactly one tier, so the two published categories
    # partition the total. A mismatch means a tier exists that neither
    # category selects, which would leave the headline count disagreeing
    # with the two beneath it on the published page.
    partitioned = counts['unit'] + counts['integration']
    if counts['total'] != partitioned:
        raise RuntimeError(
            f"badge counts do not partition: total {counts['total']} but "
            f"unit {counts['unit']} + integration {counts['integration']} "
            f'= {partitioned}. A tier marker exists that neither category '
            'selects; add it to the integration expression in _BADGES.'
        )

    for removed_name in prune_extra_files(out_dir):
        print(f'pruned: {removed_name}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
