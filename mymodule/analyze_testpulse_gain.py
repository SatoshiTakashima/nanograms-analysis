#!/usr/bin/env python3
"""Analyze NanoGRAMS VATA test-pulse data and extract time-dependent FEC gains.

The workflow is:
  1. read each configured VATA ROOT file,
  2. subtract event-by-event common-mode noise with the median over 64 channels,
  3. inspect the configured test-pulse channel for each FEC,
  4. decide whether the FEC has a usable peak,
  5. fit the usable peaks with a Gaussian,
  6. save per-FEC fit diagnostics and a compact time-series CSV.

The compact CSV is intentionally compatible with the existing interpolation
workflow: time_id,FEC0,FEC1,FEC2,FEC3.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import uproot
import yaml


@dataclass(frozen=True)
class FitConfig:
    data_root: Path
    file_name: str
    tree_name: str
    time_ids: list[str] | None
    date_dirs: list[str] | None
    test_pulse_channel: int
    fec_ids: list[int]
    fit_window_adu: float
    histogram_min_adu: float
    histogram_max_adu: float
    histogram_bin_width_adu: float
    min_events: int
    min_fit_entries: int
    min_peak_count: float
    peak_min_adu: float
    peak_max_adu: float
    sigma_min_adu: float
    sigma_max_adu: float
    max_chi2_ndf: float
    outdir: Path
    fits_csv: Path
    summary_csv: Path
    plot_dir: Path
    review_csv: Path
    review_plot_dir: Path
    make_plots: bool
    use_manual_selection: bool


@dataclass
class FitResult:
    time_id: str
    file_path: Path
    fec: int
    alive: bool = False
    peak_adu: float = math.nan
    peak_adu_err: float = math.nan
    sigma_adu: float = math.nan
    sigma_adu_err: float = math.nan
    amplitude: float = math.nan
    baseline: float = math.nan
    n_events: int = 0
    n_fit: int = 0
    chi2_ndf: float = math.nan
    fit_min_adu: float = math.nan
    fit_max_adu: float = math.nan
    reason: str = ""
    hist_counts: np.ndarray | None = None
    hist_edges: np.ndarray | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit NanoGRAMS VATA test-pulse peaks.")
    parser.add_argument("config", help="YAML configuration file.")
    return parser.parse_args()


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


def as_list(value, item_type=str) -> list | None:
    if value is None:
        return None
    if isinstance(value, str):
        if value.lower() in ("auto", "all", "none", ""):
            return None
        return [item_type(value)]
    return [item_type(item) for item in value]


def load_config(path: str | Path) -> FitConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as f:
        config = yaml.safe_load(f) or {}

    outdir = resolve_path(config_path, nested_get(config, "output.outdir", "products/testpulse_gain"))
    fits_csv = outdir / str(nested_get(config, "output.fits_csv", "testpulse_fits.csv"))
    summary_csv = outdir / str(nested_get(config, "output.summary_csv", "run10_testpulse_data.csv"))
    plot_dir = outdir / str(nested_get(config, "output.plot_dir", "plots"))
    review_csv = outdir / str(nested_get(config, "review.selection_csv", "manual_selection.csv"))
    review_plot_dir = outdir / str(nested_get(config, "review.plot_dir", "review_plots"))

    return FitConfig(
        data_root=resolve_path(config_path, required_value(config, "input.data_root")),
        file_name=str(nested_get(config, "input.file_name", "outfile00001_000.root")),
        tree_name=str(nested_get(config, "input.tree", "vatatree")),
        time_ids=as_list(nested_get(config, "input.time_ids", None), str),
        date_dirs=as_list(nested_get(config, "input.date_dirs", None), str),
        test_pulse_channel=int(nested_get(config, "test_pulse.channel", 17)),
        fec_ids=as_list(nested_get(config, "test_pulse.fec_ids", [0, 1, 2, 3]), int) or [0, 1, 2, 3],
        fit_window_adu=float(nested_get(config, "fit.fit_window_adu", 20.0)),
        histogram_min_adu=float(nested_get(config, "fit.histogram_min_adu", -50.0)),
        histogram_max_adu=float(nested_get(config, "fit.histogram_max_adu", 1024.0)),
        histogram_bin_width_adu=float(nested_get(config, "fit.histogram_bin_width_adu", 1.0)),
        min_events=int(nested_get(config, "quality.min_events", 100)),
        min_fit_entries=int(nested_get(config, "quality.min_fit_entries", 50)),
        min_peak_count=float(nested_get(config, "quality.min_peak_count", 20.0)),
        peak_min_adu=float(nested_get(config, "quality.peak_min_adu", 20.0)),
        peak_max_adu=float(nested_get(config, "quality.peak_max_adu", 1000.0)),
        sigma_min_adu=float(nested_get(config, "quality.sigma_min_adu", 0.5)),
        sigma_max_adu=float(nested_get(config, "quality.sigma_max_adu", 50.0)),
        max_chi2_ndf=float(nested_get(config, "quality.max_chi2_ndf", 5.0)),
        outdir=outdir,
        fits_csv=fits_csv,
        summary_csv=summary_csv,
        plot_dir=plot_dir,
        review_csv=review_csv,
        review_plot_dir=review_plot_dir,
        make_plots=bool(nested_get(config, "output.make_plots", True)),
        use_manual_selection=bool(nested_get(config, "review.use_manual_selection", False)),
    )


def parse_time_id(time_id: str) -> datetime:
    return datetime.strptime(time_id, "%Y%m%d/%H%M_%S")


def expand_config_time_id(time_id: str, cfg: FitConfig) -> list[str]:
    if "/" in time_id:
        return [time_id]

    date_dirs = cfg.date_dirs
    if date_dirs is None:
        date_dirs = sorted(path.name for path in cfg.data_root.iterdir() if path.is_dir())

    candidates = [f"{date_dir}/{time_id}" for date_dir in date_dirs]
    existing = [
        candidate
        for candidate in candidates
        if (cfg.data_root / candidate / cfg.file_name).exists()
    ]
    if existing:
        return existing
    if len(candidates) == 1:
        return candidates
    return candidates


def discover_time_ids(cfg: FitConfig) -> list[str]:
    if cfg.time_ids is not None:
        time_ids: list[str] = []
        for time_id in cfg.time_ids:
            time_ids.extend(expand_config_time_id(time_id, cfg))
        return sorted(time_ids, key=parse_time_id)

    date_dirs = cfg.date_dirs
    if date_dirs is None:
        date_dirs = sorted(path.name for path in cfg.data_root.iterdir() if path.is_dir())

    time_ids: list[str] = []
    for date_dir in date_dirs:
        date_path = cfg.data_root / date_dir
        if not date_path.is_dir():
            continue
        for time_path in sorted(path for path in date_path.iterdir() if path.is_dir()):
            if (time_path / cfg.file_name).exists():
                time_ids.append(f"{date_dir}/{time_path.name}")
    return sorted(time_ids, key=parse_time_id)


def gaussian_with_const(x, amplitude, mean, sigma, baseline):
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2.0) + baseline


def integer_centered_histogram_edges(cfg: FitConfig) -> np.ndarray:
    """Return bin edges whose centers sit on integer-like ADU values."""
    half_width = 0.5 * cfg.histogram_bin_width_adu
    return np.arange(cfg.histogram_min_adu - half_width,
                     cfg.histogram_max_adu + half_width + cfg.histogram_bin_width_adu,
                     cfg.histogram_bin_width_adu)


def fit_histogram_with_iminuit_mle(values: np.ndarray,
                                   cfg: FitConfig,
                                   peak_guess: float,
                                   peak_count: float,
                                   initial_mean_adu: float | None = None,
                                   initial_height: float | None = None,
                                   initial_width_adu: float | None = None,
                                   fit_window_adu: float | None = None) -> dict[str, float]:
    from iminuit import Minuit

    fit_center = (
        float(initial_mean_adu)
        if initial_mean_adu is not None and np.isfinite(initial_mean_adu)
        else peak_guess
    )
    fit_window = (
        float(fit_window_adu)
        if fit_window_adu is not None and np.isfinite(fit_window_adu) and fit_window_adu > 0.0
        else cfg.fit_window_adu
    )
    fit_min = fit_center - fit_window
    fit_max = fit_center + fit_window
    nbins = max(1, int(math.ceil((fit_max - fit_min) / cfg.histogram_bin_width_adu)))
    fit_data = values[(values >= fit_min) & (values <= fit_max)].astype(float)
    if fit_data.size == 0:
        raise RuntimeError("No entries in the fit range.")

    fit_width = fit_max - fit_min
    bin_width = fit_width / nbins
    inv_sqrt_2pi = 1.0 / math.sqrt(2.0 * math.pi)
    sqrt_2 = math.sqrt(2.0)

    def gaussian_norm(mean: float, sigma: float) -> float:
        hi = (fit_max - mean) / (sqrt_2 * sigma)
        lo = (fit_min - mean) / (sqrt_2 * sigma)
        return max(0.5 * (math.erf(hi) - math.erf(lo)), 1.0e-12)

    def pdf(x: np.ndarray, mean: float, sigma: float, signal_fraction: float) -> np.ndarray:
        norm = gaussian_norm(mean, sigma)
        z = (x - mean) / sigma
        gaussian = inv_sqrt_2pi * np.exp(-0.5 * z * z) / sigma / norm
        uniform = 1.0 / fit_width
        return signal_fraction * gaussian + (1.0 - signal_fraction) * uniform

    def nll(mean: float, sigma: float, signal_fraction: float) -> float:
        if sigma <= 0.0 or signal_fraction < 0.0 or signal_fraction > 1.0:
            return np.inf
        probabilities = pdf(fit_data, mean, sigma, signal_fraction)
        if np.any(probabilities <= 0.0) or not np.all(np.isfinite(probabilities)):
            return np.inf
        return float(-np.sum(np.log(probabilities)))

    mean0 = fit_center
    sigma0 = (
        float(initial_width_adu)
        if initial_width_adu is not None and np.isfinite(initial_width_adu) and initial_width_adu > 0.0
        else float(np.std(fit_data))
    )
    sigma0 = float(np.clip(sigma0, cfg.sigma_min_adu, cfg.sigma_max_adu))
    height0 = (
        float(initial_height)
        if initial_height is not None and np.isfinite(initial_height) and initial_height > 0.0
        else peak_count
    )
    max_gaussian_height = fit_data.size * bin_width * inv_sqrt_2pi / sigma0 / gaussian_norm(mean0, sigma0)
    signal_fraction0 = float(np.clip(height0 / max(max_gaussian_height, 1.0), 0.05, 1.0))

    minuit = Minuit(
        nll,
        mean=mean0,
        sigma=sigma0,
        signal_fraction=signal_fraction0,
    )
    minuit.errordef = Minuit.LIKELIHOOD
    minuit.limits["mean"] = (fit_min, fit_max)
    minuit.limits["sigma"] = (cfg.sigma_min_adu, cfg.sigma_max_adu)
    minuit.limits["signal_fraction"] = (0.0, 1.0)
    minuit.strategy = 1
    minuit.migrad()
    minuit.hesse()

    mean = float(minuit.values["mean"])
    sigma = abs(float(minuit.values["sigma"]))
    signal_fraction = float(np.clip(minuit.values["signal_fraction"], 0.0, 1.0))
    norm = gaussian_norm(mean, sigma)
    amplitude = fit_data.size * bin_width * signal_fraction * inv_sqrt_2pi / sigma / norm
    baseline = fit_data.size * bin_width * (1.0 - signal_fraction) / fit_width

    # The VATA ADU values are integer-like after CMN subtraction for these
    # test-pulse files. Use bins centered on integer ADU values for the chi2
    # diagnostic so it is stable against small changes of the fit-range center.
    chi_edges = np.arange(math.floor(fit_min) - 0.5,
                          math.ceil(fit_max) + 0.5 + cfg.histogram_bin_width_adu,
                          cfg.histogram_bin_width_adu)
    counts_all, chi_edges = np.histogram(values, bins=chi_edges)
    centers_all = 0.5 * (chi_edges[:-1] + chi_edges[1:])
    chi_mask = (centers_all >= fit_min) & (centers_all <= fit_max)
    counts = counts_all[chi_mask]
    centers = centers_all[chi_mask]
    expected = fit_data.size * bin_width * pdf(centers, mean, sigma, signal_fraction)
    chi2 = float(np.sum((counts - expected) ** 2 / np.maximum(expected, 1.0)))
    ndf = max(1, counts.size - 3)
    mean_err = float(minuit.errors["mean"])
    sigma_err = float(minuit.errors["sigma"])

    return {
        "valid": float(bool(minuit.valid)),
        "amplitude": float(amplitude),
        "peak_adu": mean,
        "peak_adu_err": mean_err,
        "sigma_adu": sigma,
        "sigma_adu_err": sigma_err,
        "baseline": float(baseline),
        "chi2": chi2,
        "ndf": float(ndf),
        "chi2_ndf": chi2 / ndf,
        "fit_min_adu": fit_min,
        "fit_max_adu": fit_max,
    }


def read_fec_channel_values_by_fec(file_path: Path, cfg: FitConfig) -> dict[int, np.ndarray]:
    tree = uproot.open(file_path)[cfg.tree_name]
    cut = (
        "(flag_self_trig!=0)&(good_event==1)&"
        "(hitnum==64)"
    )
    arrays = tree.arrays(["adc", "fec_index"], cut, library="np")
    values_by_fec: dict[int, np.ndarray] = {
        fec: np.empty(0, dtype=float)
        for fec in cfg.fec_ids
    }
    if len(arrays["adc"]) == 0:
        return values_by_fec

    fec_index = arrays["fec_index"]
    adc_all = arrays["adc"]
    for fec in cfg.fec_ids:
        selected_adc = adc_all[fec_index == fec]
        if len(selected_adc) == 0:
            continue
        adc = np.vstack(selected_adc).astype(float)
        cmn = np.median(adc, axis=1)
        adu_cmn = adc - cmn[:, np.newaxis]
        values_by_fec[fec] = adu_cmn[:, cfg.test_pulse_channel]
    return values_by_fec


def read_fec_channel_values(file_path: Path, cfg: FitConfig, fec: int) -> np.ndarray:
    return read_fec_channel_values_by_fec(file_path, cfg).get(fec, np.empty(0, dtype=float))


def fit_fec(time_id: str,
            file_path: Path,
            cfg: FitConfig,
            fec: int,
            initial_mean_adu: float | None = None,
            initial_height: float | None = None,
            initial_width_adu: float | None = None,
            fit_window_adu: float | None = None,
            values: np.ndarray | None = None) -> FitResult:
    if values is None:
        values = read_fec_channel_values(file_path, cfg, fec)
    result = FitResult(time_id=time_id, file_path=file_path, fec=fec, n_events=int(values.size))
    if values.size < cfg.min_events:
        result.reason = "too_few_events"
        return result

    edges = integer_centered_histogram_edges(cfg)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    result.hist_counts = counts
    result.hist_edges = edges

    peak_index = int(np.argmax(counts))
    peak_guess = float(centers[peak_index])
    peak_count = float(counts[peak_index])
    if peak_count < cfg.min_peak_count:
        result.reason = "weak_peak"
        return result

    fit_center = (
        float(initial_mean_adu)
        if initial_mean_adu is not None and np.isfinite(initial_mean_adu)
        else peak_guess
    )
    fit_window = (
        float(fit_window_adu)
        if fit_window_adu is not None and np.isfinite(fit_window_adu) and fit_window_adu > 0.0
        else cfg.fit_window_adu
    )
    fit_mask = np.abs(values - fit_center) <= fit_window
    fit_values = values[fit_mask]
    result.n_fit = int(fit_values.size)
    if fit_values.size < cfg.min_fit_entries:
        result.reason = "too_few_fit_entries"
        return result

    try:
        minuit_fit = fit_histogram_with_iminuit_mle(
            values,
            cfg,
            peak_guess,
            peak_count,
            initial_mean_adu=initial_mean_adu,
            initial_height=initial_height,
            initial_width_adu=initial_width_adu,
            fit_window_adu=fit_window,
        )
    except ModuleNotFoundError as exc:
        if exc.name == "iminuit":
            result.reason = "fit_failed:iminuit_missing"
        else:
            result.reason = f"fit_failed:{type(exc).__name__}"
        return result
    except Exception as exc:
        result.reason = f"fit_failed:{type(exc).__name__}"
        return result

    result.amplitude = float(minuit_fit["amplitude"])
    result.peak_adu = float(minuit_fit["peak_adu"])
    result.peak_adu_err = float(minuit_fit["peak_adu_err"])
    result.sigma_adu = abs(float(minuit_fit["sigma_adu"]))
    result.sigma_adu_err = float(minuit_fit["sigma_adu_err"])
    result.baseline = float(minuit_fit["baseline"])
    result.chi2_ndf = float(minuit_fit["chi2_ndf"])
    result.fit_min_adu = float(minuit_fit["fit_min_adu"])
    result.fit_max_adu = float(minuit_fit["fit_max_adu"])

    if not all(np.isfinite(value) for value in (
        result.amplitude,
        result.peak_adu,
        result.peak_adu_err,
        result.sigma_adu,
        result.sigma_adu_err,
        result.baseline,
        result.chi2_ndf,
    )):
        result.reason = "fit_invalid"
    elif not (cfg.peak_min_adu <= result.peak_adu <= cfg.peak_max_adu):
        result.reason = "peak_out_of_range"
    elif not (cfg.sigma_min_adu <= result.sigma_adu <= cfg.sigma_max_adu):
        result.reason = "sigma_out_of_range"
    elif result.chi2_ndf > cfg.max_chi2_ndf:
        result.reason = "large_chi2"
    else:
        result.alive = True
        result.reason = "ok"
    return result


def write_fits_csv(results: list[FitResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time_id", "datetime", "fec", "alive", "peak_adu", "peak_adu_err",
            "sigma_adu", "sigma_adu_err", "amplitude", "baseline", "n_events",
            "n_fit", "chi2_ndf", "fit_min_adu", "fit_max_adu", "reason", "file_path",
        ])
        for result in results:
            writer.writerow([
                result.time_id,
                parse_time_id(result.time_id).isoformat(sep=" "),
                result.fec,
                int(result.alive),
                result.peak_adu,
                result.peak_adu_err,
                result.sigma_adu,
                result.sigma_adu_err,
                result.amplitude,
                result.baseline,
                result.n_events,
                result.n_fit,
                result.chi2_ndf,
                result.fit_min_adu,
                result.fit_max_adu,
                result.reason,
                result.file_path,
            ])


def parse_manual_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text == "":
        return None
    if text in ("1", "true", "t", "yes", "y", "ok"):
        return True
    if text in ("0", "false", "f", "no", "n", "ng"):
        return False
    raise ValueError(f"manual selection must be 1/0, true/false, ok/ng, or blank: {value}")


def load_manual_selection(path: Path) -> dict[tuple[str, int], bool]:
    if not path.exists():
        return {}

    selection: dict[tuple[str, int], bool] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            use_fit = parse_manual_bool(row.get("use_fit"))
            if use_fit is None:
                continue
            selection[(row["time_id"], int(row["fec"]))] = use_fit
    return selection


def result_is_selected(result: FitResult,
                       cfg: FitConfig,
                       manual_selection: dict[tuple[str, int], bool]) -> bool:
    if not cfg.use_manual_selection:
        return result.alive
    return bool(manual_selection.get((result.time_id, result.fec), False))


def missing_manual_selection(results: list[FitResult],
                             manual_selection: dict[tuple[str, int], bool]) -> list[tuple[str, int]]:
    return [
        (result.time_id, result.fec)
        for result in results
        if (result.time_id, result.fec) not in manual_selection
    ]


def write_manual_review_csv(results: list[FitResult],
                            cfg: FitConfig,
                            manual_selection: dict[tuple[str, int], bool]) -> None:
    cfg.review_csv.parent.mkdir(parents=True, exist_ok=True)
    with cfg.review_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time_id", "datetime", "fec", "use_fit", "suggested_use_fit",
            "peak_adu", "peak_adu_err", "sigma_adu", "sigma_adu_err", "chi2_ndf",
            "n_events", "fit_min_adu", "fit_max_adu", "reason", "review_plot",
            "fit_plot", "file_path",
        ])
        for result in results:
            key = (result.time_id, result.fec)
            use_fit = "" if key not in manual_selection else int(manual_selection[key])
            safe_time = result.time_id.replace("/", "_")
            writer.writerow([
                result.time_id,
                parse_time_id(result.time_id).isoformat(sep=" "),
                result.fec,
                use_fit,
                int(result.alive),
                result.peak_adu,
                result.peak_adu_err,
                result.sigma_adu,
                result.sigma_adu_err,
                result.chi2_ndf,
                result.n_events,
                result.fit_min_adu,
                result.fit_max_adu,
                result.reason,
                cfg.review_plot_dir / f"{safe_time}.png",
                cfg.plot_dir / f"{safe_time}_FEC{result.fec}.png",
                result.file_path,
            ])


def write_summary_csv(results: list[FitResult],
                      cfg: FitConfig,
                      manual_selection: dict[tuple[str, int], bool]) -> None:
    by_time: dict[str, dict[int, FitResult]] = {}
    for result in results:
        by_time.setdefault(result.time_id, {})[result.fec] = result

    cfg.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with cfg.summary_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_id", "datetime", *[f"FEC{fec}" for fec in cfg.fec_ids]])
        for time_id in sorted(by_time, key=parse_time_id):
            row = [time_id, parse_time_id(time_id).isoformat(sep=" ")]
            for fec in cfg.fec_ids:
                result = by_time[time_id].get(fec)
                row.append(
                    result.peak_adu
                    if result is not None and result_is_selected(result, cfg, manual_selection)
                    else float("nan")
                )
            writer.writerow(row)


def plot_fit(result: FitResult, cfg: FitConfig) -> None:
    if result.hist_counts is None or result.hist_edges is None:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg.plot_dir.mkdir(parents=True, exist_ok=True)
    centers = 0.5 * (result.hist_edges[:-1] + result.hist_edges[1:])
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.stairs(result.hist_counts, result.hist_edges, label=f"FEC{result.fec}")
    if np.isfinite(result.peak_adu) and np.isfinite(result.sigma_adu):
        fit_min = result.fit_min_adu if np.isfinite(result.fit_min_adu) else result.peak_adu - cfg.fit_window_adu
        fit_max = result.fit_max_adu if np.isfinite(result.fit_max_adu) else result.peak_adu + cfg.fit_window_adu
        x = np.linspace(fit_min, fit_max, 400)
        ax.plot(x,
                gaussian_with_const(x, result.amplitude, result.peak_adu,
                                    result.sigma_adu, result.baseline),
                color="tab:red",
                label=f"fit: {result.peak_adu:.2f} ADU")
    ax.set_title(f"{result.time_id} FEC{result.fec}: {result.reason}")
    ax.set_xlabel(f"ch{cfg.test_pulse_channel} ADU - CMN")
    ax.set_ylabel("Counts")
    ax.set_xlim(cfg.histogram_min_adu, cfg.histogram_max_adu)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    safe_time = result.time_id.replace("/", "_")
    fig.savefig(cfg.plot_dir / f"{safe_time}_FEC{result.fec}.png", dpi=140)
    plt.close(fig)


def plot_review_page(time_id: str, results: list[FitResult], cfg: FitConfig) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg.review_plot_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
    axes_flat = axes.ravel()

    for ax, result in zip(axes_flat, sorted(results, key=lambda item: item.fec)):
        if result.hist_counts is not None and result.hist_edges is not None:
            centers = 0.5 * (result.hist_edges[:-1] + result.hist_edges[1:])
            ax.stairs(result.hist_counts, result.hist_edges, label=f"FEC{result.fec}")
            if np.isfinite(result.peak_adu) and np.isfinite(result.sigma_adu):
                fit_min = result.fit_min_adu if np.isfinite(result.fit_min_adu) else result.peak_adu - cfg.fit_window_adu
                fit_max = result.fit_max_adu if np.isfinite(result.fit_max_adu) else result.peak_adu + cfg.fit_window_adu
                x = np.linspace(fit_min, fit_max, 400)
                ax.plot(x,
                        gaussian_with_const(x, result.amplitude, result.peak_adu,
                                            result.sigma_adu, result.baseline),
                        color="tab:red",
                        lw=1.8,
                        label=f"{result.peak_adu:.2f} ADU")
        else:
            ax.text(0.5, 0.5, "No histogram", ha="center", va="center",
                    transform=ax.transAxes)

        suggestion = "suggest OK" if result.alive else "suggest NG"
        ax.set_title(f"FEC{result.fec}: {suggestion}, {result.reason}")
        ax.set_xlim(cfg.histogram_min_adu, cfg.histogram_max_adu)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

    for ax in axes[-1, :]:
        ax.set_xlabel(f"ch{cfg.test_pulse_channel} ADU - CMN")
    for ax in axes[:, 0]:
        ax.set_ylabel("Counts")

    fig.suptitle(f"{time_id} test-pulse review")
    fig.tight_layout()
    safe_time = time_id.replace("/", "_")
    fig.savefig(cfg.review_plot_dir / f"{safe_time}.png", dpi=150)
    plt.close(fig)


def plot_time_trend(results: list[FitResult],
                    cfg: FitConfig,
                    manual_selection: dict[tuple[str, int], bool]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfg.plot_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for fec in cfg.fec_ids:
        fec_results = [
            result for result in results
            if result.fec == fec and result_is_selected(result, cfg, manual_selection)
        ]
        if not fec_results:
            continue
        x = [parse_time_id(r.time_id) for r in fec_results]
        y = [r.peak_adu for r in fec_results]
        yerr = [r.peak_adu_err for r in fec_results]
        ax.errorbar(x, y, yerr=yerr, marker="o", ls="-", capsize=2, label=f"FEC{fec}")
    ax.set_xlabel("Time")
    ax.set_ylabel(f"ch{cfg.test_pulse_channel} peak (ADU - CMN)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(cfg.plot_dir / "time_trend.png", dpi=160)
    plt.close(fig)


def main() -> None:
    cfg = load_config(parse_args().config)
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    if cfg.make_plots:
        mpl_cache_dir = cfg.outdir / ".matplotlib"
        mpl_cache_dir.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache_dir))
        cfg.plot_dir.mkdir(parents=True, exist_ok=True)
        cfg.review_plot_dir.mkdir(parents=True, exist_ok=True)

    results: list[FitResult] = []
    for time_id in discover_time_ids(cfg):
        file_path = cfg.data_root / time_id / cfg.file_name
        if not file_path.exists():
            print(f"[skip] missing: {file_path}")
            continue
        print(f"[read] {time_id}: {file_path}")
        values_by_fec = read_fec_channel_values_by_fec(file_path, cfg)
        for fec in cfg.fec_ids:
            result = fit_fec(time_id, file_path, cfg, fec, values=values_by_fec.get(fec))
            results.append(result)
            print(
                f"  FEC{fec}: {result.reason}, "
                f"peak={result.peak_adu:.3g}, sigma={result.sigma_adu:.3g}, "
                f"chi2/ndf={result.chi2_ndf:.3g}, n={result.n_events}"
            )
            if cfg.make_plots:
                plot_fit(result, cfg)

    if cfg.make_plots:
        by_time: dict[str, list[FitResult]] = {}
        for result in results:
            by_time.setdefault(result.time_id, []).append(result)
        for time_id, time_results in by_time.items():
            plot_review_page(time_id, time_results, cfg)

    manual_selection = load_manual_selection(cfg.review_csv)
    write_fits_csv(results, cfg.fits_csv)
    write_manual_review_csv(results, cfg, manual_selection)

    missing_selection = missing_manual_selection(results, manual_selection)
    if cfg.use_manual_selection and missing_selection:
        print(
            f"[review] Edit use_fit in {cfg.review_csv}, then rerun this script. "
            f"{len(missing_selection)} FEC entries are still blank."
        )
    else:
        write_summary_csv(results, cfg, manual_selection)
        if cfg.make_plots:
            plot_time_trend(results, cfg, manual_selection)
        print(f"summary: {cfg.summary_csv}")

    print(f"fits:    {cfg.fits_csv}")
    print(f"review:  {cfg.review_csv}")
    if cfg.make_plots:
        print(f"plots:   {cfg.plot_dir}")
        print(f"review plots: {cfg.review_plot_dir}")


if __name__ == "__main__":
    main()
