#!/usr/bin/env python3
"""Interpolate test-pulse FEC gains for gamma-ray run times.

Input test-pulse CSV is expected to contain:
  time_id,datetime,FEC0,FEC1,FEC2,FEC3

The output CSV keeps the gamma run identifiers and appends interpolated FEC
peak ADUs.  These values can be passed to NanoGRAMSCalibration as gain_tp_hash.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

def nested_get(config: dict[str, Any], path: str, default=None):
    node: Any = config
    for key in path.split("."):
        if not isinstance(node, dict) or key not in node:
            return default
        node = node[key]
    return node


def required_value(config: dict[str, Any], path: str):
    value = nested_get(config, path)
    if value is None:
        raise KeyError(f"Missing required config value: {path}")
    return value


def resolve_path(config_path: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def load_config(path: str | Path) -> tuple[Path, Path, Path, str, list[int]]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as f:
        config = yaml.safe_load(f) or {}

    testpulse_csv = resolve_path(config_path, required_value(config, "input.testpulse_csv"))
    gamma_csv = resolve_path(config_path, required_value(config, "input.gamma_csv"))
    output_csv = resolve_path(config_path, nested_get(config, "output.csv", "gamma_runs_with_gain.csv"))
    gamma_time_column = str(nested_get(config, "input.gamma_time_column", "time"))
    fec_ids = [int(fec) for fec in nested_get(config, "input.fec_ids", [0, 1, 2, 3])]
    return testpulse_csv, gamma_csv, output_csv, gamma_time_column, fec_ids


def normalize_gamma_time_ids(series: pd.Series) -> pd.Series:
    values = series.astype(str)
    missing_date = ~values.str.contains("/")
    if missing_date.any():
        examples = ", ".join(values[missing_date].head(3))
        raise ValueError(
            "Gamma time values must be YYYYMMDD/HHMM_SS. "
            f"Missing date in: {examples}"
        )
    return values


def parse_gamma_datetime(series: pd.Series) -> pd.Series:
    return pd.to_datetime(normalize_gamma_time_ids(series), format="%Y%m%d/%H%M_%S")


def interpolate_column(tp: pd.DataFrame, gamma_times: pd.Series, column: str) -> np.ndarray:
    """Interpolate one FEC gain column at gamma-ray run times.

    Each FEC is handled independently.  A gamma run before the first valid
    test-pulse gain for that FEC is left blank instead of being filled from a
    future measurement.  Values are written only when the gamma time is exactly
    on a valid test-pulse point or bracketed by valid points on both sides.
    """
    valid = tp[["datetime", column]].dropna().sort_values("datetime")
    if valid.empty:
        return np.full(len(gamma_times), np.nan)

    x = valid["datetime"].astype("int64").to_numpy(dtype=float) / 1.0e9
    y = valid[column].to_numpy(dtype=float)
    x_gamma = gamma_times.astype("int64").to_numpy(dtype=float) / 1.0e9

    result = np.full(len(x_gamma), np.nan, dtype=float)
    index_after = np.searchsorted(x, x_gamma, side="left")

    exact = np.zeros(len(x_gamma), dtype=bool)
    has_after = index_after < len(x)
    exact[has_after] = x[index_after[has_after]] == x_gamma[has_after]
    result[exact] = y[index_after[exact]]

    bracketed = (~exact) & (index_after > 0) & (index_after < len(x))
    before = index_after[bracketed] - 1
    after = index_after[bracketed]
    weight_after = (x_gamma[bracketed] - x[before]) / (x[after] - x[before])
    result[bracketed] = (1.0 - weight_after) * y[before] + weight_after * y[after]
    return result
