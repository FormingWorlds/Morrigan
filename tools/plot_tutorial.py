"""Draw the figure the standalone tutorial shows, in both colour schemes.

Reads a finished run's ``full_systems`` table and plots each body's semi-major
axis and mass against time, so the mergers are visible as tracks that end and a
survivor that steps up in mass.

Styling comes from proteus-mpl, the PROTEUS matplotlib theme. It is used through
the installed package when that is available, and read from a checkout of the
visual-language repository otherwise, so the figure looks the same either way.

Each figure carries the version of the model that produced the run it plots, so
a reader can tell what the figure is a figure of. That is the installed
``morrigan`` version, which for a released install is the release tag itself,
not the commit of whatever checkout happened to draw the figure.

Usage: python tools/plot_tutorial.py <run_directory> <output_directory> [--label TEXT]
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

M_EARTH = 5.972e24
YEAR = 365 * 24 * 3600

# Fallback location of the theme when the package is not installed.
VL_STYLE = Path.home() / 'git/proteus-visual-language/figures/proteus-mpl/src/proteus_mpl'

# proteus-mpl's own cycles, neutral first.
CYCLE = ['#10151B', '#E23D28', '#E0A32E', '#57A05C', '#1B6FA8', '#593E74']
CYCLE_DARK = ['#F2F5F7', '#E23D28', '#E0A32E', '#57C08A', '#A8D4E8', '#9B7BBE']

# Ink for the provenance label, secondary text in each scheme.
LABEL_INK = {'light': '#7A8894', 'dark': '#5A6B7A'}


def use_theme(variant: str) -> None:
    """Apply proteus-mpl, from the package if installed, else from a checkout."""
    try:
        import proteus_mpl

        proteus_mpl.use(variant)
    except ImportError:
        base = VL_STYLE / 'proteus.mplstyle'
        if not base.exists():
            raise SystemExit(
                'proteus-mpl is not installed and no visual-language checkout was '
                f'found at {VL_STYLE}. Install it with: pip install proteus-mpl'
            ) from None
        if variant == 'dark':
            plt.style.use([str(base), str(VL_STYLE / 'proteus_dark.mplstyle')])
            mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CYCLE_DARK)
        else:
            plt.style.use(str(base))
            mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CYCLE)
    # the documentation places these on the page background of either scheme
    for key in ('figure.facecolor', 'savefig.facecolor', 'axes.facecolor'):
        mpl.rcParams[key] = 'none'
    mpl.rcParams['savefig.transparent'] = True


def model_label() -> str:
    """Identify the model that produced the run, not the checkout drawing it."""
    try:
        from morrigan import __version__ as version
    except Exception:
        return 'morrigan (version unknown)'
    if version in ('0.0.0', ''):
        # setuptools-scm fallback: a checkout with no tag history in view
        return 'morrigan (untagged checkout)'
    dirty = '+dirty' if version.endswith('.dirty') or '.d2' in version else ''
    sha = tag_commit(version)
    return f'morrigan {version}{dirty}' + (f' ({sha})' if sha else '')


def tag_commit(version: str) -> str:
    """The commit a released version points at, when this is a checkout of it.

    setuptools-scm reports the tag 26.07.26 as the version 26.7.26, so the
    zero-padded tag has to be rebuilt before it can be looked up.
    """
    parts = version.split('.')
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        return ''
    tag = '.'.join(f'{int(part):02d}' for part in parts)
    try:
        out = subprocess.run(
            ['git', 'rev-list', '-n1', '--abbrev-commit', tag],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return ''
    return out.stdout.strip() if out.returncode == 0 else ''


def draw(run: Path, out: Path, variant: str, label: str, system: int = 0) -> Path:
    use_theme(variant)
    cycle = CYCLE_DARK if variant == 'dark' else CYCLE
    source = run / f'data/full_systems/full_system_{system:02d}.csv'
    if not source.exists():
        available = sorted(p.name for p in (run / 'data/full_systems').glob('*.csv'))
        raise SystemExit(
            f'no such system table: {source}. '
            + (f'This run has {available}.' if available else 'This run has none.')
        )
    table = ascii.read(source, format='fixed_width')
    bodies = sorted(set(table['id']))
    t_end = np.asarray(table['t']).max() / YEAR

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'hspace': 0.14})

    for n, body in enumerate(bodies):
        track = table[table['id'] == body]
        t_yr = np.asarray(track['t']) / YEAR
        colour = cycle[n % len(cycle)]
        axes[0].plot(t_yr, track['a_AU'], color=colour, lw=2.0, label=f'{body}')
        axes[1].plot(t_yr, np.asarray(track['Mp']) / M_EARTH, color=colour, lw=2.0)
        if len(t_yr) and t_yr[-1] < t_end:  # merged away: mark where the track ends
            axes[0].plot(t_yr[-1], track['a_AU'][-1], 'o', ms=6, color=colour)
            axes[1].plot(t_yr[-1], track['Mp'][-1] / M_EARTH, 'o', ms=6, color=colour)

    axes[0].set_ylabel('semi-major axis [AU]')
    axes[1].set_ylabel(r'mass [M$_\oplus$]')
    axes[1].set_xlabel('time [yr]')
    axes[0].set_xscale('log')

    # The legend sits outside the axes, so no track can run underneath it.
    axes[0].legend(
        title='body',
        frameon=False,
        ncol=len(bodies),
        fontsize=11,
        title_fontsize=11,
        loc='lower center',
        bbox_to_anchor=(0.5, 1.02),
        borderaxespad=0.0,
        columnspacing=1.4,
        handlelength=1.6,
    )

    # Headroom above each track, so nothing touches the frame or the legend.
    for ax in axes:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.12 * (hi - lo))

    fig.text(
        0.995,
        0.005,
        label,
        ha='right',
        va='bottom',
        fontsize=8,
        color=LABEL_INK[variant],
        family='monospace',
    )

    out.mkdir(parents=True, exist_ok=True)
    path = out / f'tutorial_standalone_{variant}.png'
    fig.savefig(path, dpi=200, bbox_inches='tight', transparent=True)
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run')
    ap.add_argument('out')
    ap.add_argument(
        '--label', default=None, help='provenance text for the corner of the figure'
    )
    ap.add_argument('--system', type=int, default=0, help='which system of a batch to plot')
    args = ap.parse_args()
    label = args.label or model_label()
    for variant in ('light', 'dark'):
        print('wrote', draw(Path(args.run), Path(args.out), variant, label, args.system))


if __name__ == '__main__':
    main()
