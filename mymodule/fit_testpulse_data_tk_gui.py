#!/usr/bin/env python3
"""Compatibility wrapper for the customtkinter test-pulse review GUI."""

try:
    from .fit_testpulse_data_gui import *  # noqa: F401,F403
except ImportError:
    from fit_testpulse_data_gui import *  # noqa: F401,F403
