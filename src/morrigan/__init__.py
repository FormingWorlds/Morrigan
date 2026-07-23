"""
!!! info "`morrigan`"
    Semi-analytical model for the dynamical evolution of planetary systems via giant impacts, following Kimura et al. (2025).
    Author(s): Anna Grace Ulses
"""

from __future__ import annotations

try:
    from ._version import __version__, __version_tuple__
except ImportError:
    # Fallback for when the package is not installed (e.g., running from
    # source without setuptools-scm having generated _version.py).
    __version__ = '0.0.0.dev0'
    __version_tuple__ = (0, 0, 0, 'dev0')

__all__ = ['run_system', '__version__', '__version_tuple__']


def __getattr__(name):
    #import the model lazily so `import morrigan` (e.g. to read __version__) does
    #not pull in astropy and the rest of the model; a caller that reaches for
    #run_system, including hasattr(morrigan, 'run_system'), triggers this
    if name == 'run_system':
        from morrigan.driver import run_system
        return run_system
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
