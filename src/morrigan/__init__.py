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
