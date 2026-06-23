"""Compatibility imports for old notebooks.

Backprojection utilities now live in ``backprojection.py``.
"""

try:
    from .backprojection import *  # noqa: F401,F403
except ImportError:
    from backprojection import *  # noqa: F401,F403
