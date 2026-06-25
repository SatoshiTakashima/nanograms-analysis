#!/usr/bin/env python3
"""Interpolate test-pulse FEC gains for gamma-ray run times.

Input test-pulse CSV is expected to contain:
  time_id,datetime,FEC0,FEC1,FEC2,FEC3

The output CSV keeps the gamma run identifiers and appends interpolated FEC
peak ADUs.  These values can be passed to NanoGRAMSCalibration as gain_tp_hash.
"""

from __future__ import annotations

import pandas as pd
import sys
sys.path.append("../../mymodule")
from gain_interpolation import (
    load_config,
    normalize_gamma_time_ids,
    parse_gamma_datetime,
    interpolate_column,
)


if __name__ == "__main__":
    testpulse_csv, gamma_csv, output_csv, gamma_time_column, fec_ids = load_config(
        "../metadata/config_interpolate_gain.yaml"
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

    gamma[gamma_time_column] = normalize_gamma_time_ids(gamma[gamma_time_column])
    gamma["datetime"] = parse_gamma_datetime(gamma[gamma_time_column])

    for fec in fec_ids:
        column = f"FEC{fec}"
        gamma[column] = interpolate_column(tp, gamma["datetime"], column)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    gamma.to_csv(output_csv, index=False)
    print(f"saved: {output_csv}")
