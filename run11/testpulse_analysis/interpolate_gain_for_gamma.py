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


#def parse_args() -> argparse.Namespace:
#    parser = argparse.ArgumentParser(description="Interpolate test-pulse gains for gamma runs.")
#    parser.add_argument("config", help="YAML configuration file.")
#    return parser.parse_args()


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


def load_config(path: str | Path) -> tuple[Path, Path, Path, str, str | None, list[int]]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as f:
        config = yaml.safe_load(f) or {}

    testpulse_csv = resolve_path(config_path, required_value(config, "input.testpulse_csv"))
    gamma_csv = resolve_path(config_path, required_value(config, "input.gamma_csv"))
    output_csv = resolve_path(config_path, nested_get(config, "output.csv", "gamma_runs_with_gain.csv"))
    gamma_time_column = str(nested_get(config, "input.gamma_time_column", "time"))
    default_date = nested_get(config, "input.default_date", None)
    fec_ids = [int(fec) for fec in nested_get(config, "input.fec_ids", [0, 1, 2, 3])]
    return testpulse_csv, gamma_csv, output_csv, gamma_time_column, default_date, fec_ids


def parse_gamma_datetime(series: pd.Series, default_date: str | None) -> pd.Series:
    values = series.astype(str)
    if values.str.contains("/").any():
        return pd.to_datetime(values, format="%Y%m%d/%H%M_%S")
    if default_date is None:
        raise ValueError("gamma times have no date; set input.default_date in YAML.")
    return pd.to_datetime(default_date + "/" + values, format="%Y%m%d/%H%M_%S")


def interpolate_column(tp: pd.DataFrame, gamma_times: pd.Series, column: str) -> np.ndarray:
    valid = tp[["datetime", column]].dropna().sort_values("datetime")
    if valid.empty:
        return np.full(len(gamma_times), np.nan)

    x = valid["datetime"].astype("int64").to_numpy(dtype=float) / 1.0e9
    y = valid[column].to_numpy(dtype=float)
    x_gamma = gamma_times.astype("int64").to_numpy(dtype=float) / 1.0e9

    result = np.interp(x_gamma, x, y)
    result[x_gamma < x[0]] = np.nan
    result[x_gamma > x[-1]] = np.nan
    return result


    #def main() -> None:
    #    testpulse_csv, gamma_csv, output_csv, gamma_time_column, default_date, fec_ids = load_config(
    #        parse_args().config
    #    )
    #
    #    tp = pd.read_csv(testpulse_csv)
    #    if "datetime" in tp.columns:
    #        tp["datetime"] = pd.to_datetime(tp["datetime"])
    #    else:
    #        tp["datetime"] = pd.to_datetime(tp["time_id"], format="%Y%m%d/%H%M_%S")
    #    tp = tp.sort_values("datetime").reset_index(drop=True)
    #
    #    gamma = pd.read_csv(gamma_csv)
    #    if gamma_time_column not in gamma.columns:
    #        raise KeyError(f"Missing gamma time column: {gamma_time_column}")
    #    gamma["datetime"] = parse_gamma_datetime(gamma[gamma_time_column], default_date)
    #
    #    for fec in fec_ids:
    #        column = f"FEC{fec}"
    #        gamma[column] = interpolate_column(tp, gamma["datetime"], column)
    #
    #    output_csv.parent.mkdir(parents=True, exist_ok=True)
    #    gamma.to_csv(output_csv, index=False)
    #    print(f"saved: {output_csv}")


if __name__ == "__main__":
    #main()
    testpulse_csv, gamma_csv, output_csv, gamma_time_column, default_date, fec_ids = load_config(
        "../metadata/config_interpolate_gain.yaml"
        #parse_args().config
    )

    tp = pd.read_csv(testpulse_csv)
    if "datetime" in tp.columns:
        tp["datetime"] = pd.to_datetime(tp["datetime"])
    else:
        tp["datetime"] = pd.to_datetime(tp["time_id"], format="%Y%m%d/%H%M_%S")
    tp = tp.sort_values("datetime").reset_index(drop=True)

    gamma = pd.read_csv(gamma_csv)
    if gamma_time_column not in gamma.columns:
        raise KeyError(f"Missing gamma time column: {gamma_time_column}")
    gamma["datetime"] = parse_gamma_datetime(gamma[gamma_time_column], default_date)

    for fec in fec_ids:
        column = f"FEC{fec}"
        gamma[column] = interpolate_column(tp, gamma["datetime"], column)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    gamma.to_csv(output_csv, index=False)
    print(f"saved: {output_csv}")
