#!/usr/bin/env python3
import pandas as pd
from datetime import datetime

# 入力ファイル
testpulse_file = "run10_testpulse_data.csv"
data_file = "run10data.csv"

# 出力ファイル
output_file = "run10data_with_interpolated_FEC.csv"

# -----------------------------
# run10_testpulse_data.csv を読む
# -----------------------------
df_tp = pd.read_csv(testpulse_file)

# time_id: 20251001/0124_04 形式を datetime に変換
df_tp["datetime"] = pd.to_datetime(
    df_tp["time_id"],
    format="%Y%m%d/%H%M_%S"
)

# 念のため時刻順に並べる
df_tp = df_tp.sort_values("datetime").reset_index(drop=True)

# -----------------------------
# run10data.csv を読む
# -----------------------------
df_data = pd.read_csv(data_file, header=None, names=["time"])

# run10data.csv には日付がないので、
# run10_testpulse_data.csv の最初の日付を使う
date_str = df_tp["datetime"].dt.strftime("%Y%m%d").iloc[0]

# time: 0128_54 形式を datetime に変換
df_data["datetime"] = pd.to_datetime(
    date_str + "/" + df_data["time"],
    format="%Y%m%d/%H%M_%S"
)

# -----------------------------
# 時間方向の線形補間
# -----------------------------
fec_cols = ["FEC0", "FEC1", "FEC2", "FEC3"]

for fec in fec_cols:
    interpolated_values = []

    for t in df_data["datetime"]:
        # t より前の testpulse
        before = df_tp[df_tp["datetime"] <= t].tail(1)

        # t より後の testpulse
        after = df_tp[df_tp["datetime"] >= t].head(1)

        if before.empty or after.empty:
            # 範囲外なら NaN にする
            interpolated_values.append(float("nan"))
            continue

        t0 = before["datetime"].iloc[0]
        t1 = after["datetime"].iloc[0]

        y0 = before[fec].iloc[0]
        y1 = after[fec].iloc[0]

        if t0 == t1:
            # ちょうど testpulse の時刻と一致した場合
            y = y0
        else:
            # 時間差に基づく重み
            w = (t - t0).total_seconds() / (t1 - t0).total_seconds()

            # 線形補間
            y = (1 - w) * y0 + w * y1

        interpolated_values.append(y)

    df_data[fec] = interpolated_values

# 不要なら datetime 列は消す
df_out = df_data.drop(columns=["datetime"])

# CSV に保存
df_out.to_csv(output_file, index=False)

print(f"Saved: {output_file}")
