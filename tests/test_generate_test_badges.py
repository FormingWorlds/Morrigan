"""Tests for ``tools/generate_test_badges.py``, the badge count generator.

The tool publishes numbers a public page renders, and it reads them out of
pytest's human-readable summary line, so the two things that can silently
corrupt a badge are the summary regex and the exit-code handling. Both are
pinned here against the real shapes pytest emits, with the subprocess call
stubbed so no suite is collected: the point is to fix the parsing contract,
not to re-measure the suite.

See also:
- docs/How-to/testing.md
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))

import generate_test_badges as gtb  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _completed(stdout: str, returncode: int = 0) -> subprocess.CompletedProcess:
    """A stand-in for the pytest collection subprocess."""
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr='')


@pytest.mark.parametrize(
    ('summary', 'expected', 'wrong'),
    [
        # Unfiltered collection: pytest reports a bare total.
        ('39 tests collected in 0.75s\n', 39, (75,)),
        # Marker-filtered: the selected count leads, deselected trails. The
        # tool must take the numerator, not the total and not the deselected.
        ('26/39 tests collected (13 deselected) in 0.80s\n', 26, (39, 13, 80)),
        # A single test is reported in the singular.
        ('1 test collected in 0.10s\n', 1, (10,)),
    ],
    ids=['unfiltered', 'marker-filtered', 'singular'],
)
def test_the_selected_count_is_read_from_each_summary_shape(
    monkeypatch, summary, expected, wrong
):
    """The count comes from the selected figure in every summary pytest emits.

    The marker-filtered case is the discriminating one: it carries four
    numbers, and taking the wrong one yields 39, 13 or 80 rather than 26,
    every one of them a plausible-looking badge value. Each case therefore
    also asserts the result is none of the other numbers on the line.
    """
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: _completed(summary))
    count = gtb.count_tests('unit and not skip')

    assert count == expected
    # Discrimination: none of the other integers present in the summary.
    for other in wrong:
        assert count != other


def test_an_empty_selection_is_zero_rather_than_a_failure(monkeypatch, tmp_path):
    """Exit code 5 means nothing matched, which is a real count of zero.

    Morrigan has no ``slow`` tests, so this path runs on every real
    invocation; treating it as an error would break the tool outright.
    """
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: _completed('', returncode=5))
    count = gtb.count_tests('slow and not skip')

    assert count == 0
    # The zero must survive into the badge as the string "0"; shields renders
    # a missing or non-string message as a broken image rather than a zero.
    payload = json.loads(gtb.write_badge(tmp_path, 'slow', 'slow tests', count).read_text())
    assert payload['message'] == '0'


def test_a_broken_collection_raises_instead_of_publishing_a_short_count(monkeypatch):
    """A test module that fails to import must abort, not publish silently.

    This is the single most valuable behaviour in the tool: pytest still
    prints a summary when a module fails to import, so parsing on regardless
    would publish a count short by however many tests that module held.
    """
    broken = '35 tests collected, 1 error in 0.90s\n'
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: _completed(broken, returncode=2))
    with pytest.raises(RuntimeError, match='exited with'):
        gtb.count_tests('not skip')

    # An unparseable summary is also a hard failure rather than a guess.
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: _completed('no summary here\n'))
    with pytest.raises(RuntimeError, match='summary line not found'):
        gtb.count_tests('not skip')


def test_the_written_badge_is_a_valid_shields_endpoint(tmp_path):
    """Each file carries the four keys shields requires, with the count as text.

    ``message`` must be a string: shields rejects a bare integer, and the
    failure surfaces only as a broken image on the published page.
    """
    written = gtb.write_badge(tmp_path, 'unit', 'unit tests', 26)
    payload = json.loads(written.read_text())

    assert written.name == 'tests-unit.json'
    assert payload['schemaVersion'] == 1
    assert payload['label'] == 'unit tests'
    assert payload['message'] == '26'
    assert isinstance(payload['message'], str)
    assert payload['color'] == 'blue'
