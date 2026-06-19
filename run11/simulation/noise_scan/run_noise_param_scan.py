#!/usr/bin/env python3
"""Scan detector noise_level param1 and compare simulation hittrees to data.

For each configured param1 value this script:
  1. copies detector_parameters.xml and edits <noise_level param1="...">,
  2. optionally runs the configured Ruby simulation script with that XML,
  3. reads the produced simulation hittree,
  4. compares the 1-hit energy histogram to the experimental hittree,
  5. saves reduced chi-square results and best-fit plots.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import uproot
import yaml


@dataclass(frozen=True)
class ScanConfig:
    outdir: Path
    execute: bool
    num_events: int
    random_seed: int
    energy_kev: float
    ruby_command: str
    simulation_runner: Path
    detector_parameters_template: Path
    param1_values: list[float]
    sim_output_name: str
    sim_hittree_name: str
    postprocess_command: str | None
    data_hittree: Path
    tree_name: str
    energy_branch: str
    event_branch: str
    histogram_min: float
    histogram_max: float
    bin_width: float
    compare_min: float
    compare_max: float
    normalize_min: float
    spectrum_mode: str
    one_hit_only: bool


@dataclass(frozen=True)
class ScanResult:
    param1: float
    sim_hittree: Path
    alpha: float
    chi2: float
    ndf: int
    n_data: int
    n_sim: int

    @property
    def chi2_ndf(self) -> float:
        return self.chi2 / self.ndf if self.ndf > 0 else float("nan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan detector noise param1.")
    parser.add_argument("config", help="YAML configuration file.")
    return parser.parse_args()


def resolve_path(config_path: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def load_config(path: str | Path) -> ScanConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as f:
        config = yaml.safe_load(f)
    simulation = config["simulation"]
    scan = config["scan"]
    data = config["data"]
    branches = config["branches"]
    spectrum = config["spectrum"]
    histogram_node = config["histogram"]
    comparison = config["comparison"]
    output = config["output"]
    param1_values = (
        [float(value) for value in scan["values"]]
        if "values" in scan
        else np.linspace(float(scan["min"]), float(scan["max"]), int(scan["steps"])).tolist()
    )

    return ScanConfig(
        outdir              = resolve_path(config_path, output["outdir"]),
        execute             = bool(simulation["execute"]),
        num_events          = int(simulation["num_events"]),
        random_seed         = int(simulation["random_seed"]),
        energy_kev          = float(simulation["energy_kev"]),
        ruby_command        = str(simulation["ruby_command"]),
        simulation_runner   = resolve_path(config_path, simulation["runner"]),
        detector_parameters_template=resolve_path(config_path, simulation["detector_parameters_template"]),
        param1_values       = param1_values,
        sim_output_name     = str(simulation["output_root"]),
        sim_hittree_name    = str(simulation["hittree_root"]),
        postprocess_command = simulation["postprocess_command"],
        data_hittree        = resolve_path(config_path, data["hittree"]),
        tree_name           = str(data["tree"]),
        energy_branch       = str(branches["energy"]),
        event_branch        = str(branches["event"]),
        histogram_min       = float(histogram_node["min_kev"]),
        histogram_max       = float(histogram_node["max_kev"]),
        bin_width           = float(histogram_node["bin_width_kev"]),
        compare_min         = float(comparison["min_kev"]),
        compare_max         = float(comparison["max_kev"]),
        normalize_min       = float(comparison["normalize_min_kev"]),
        spectrum_mode       = str(spectrum["mode"]),
        one_hit_only        = bool(spectrum["one_hit_only"]),
    )


def run_label(param1: float) -> str:
    return f"param1_{param1:.6g}".replace(".", "p").replace("-", "m")


def write_detector_parameters(template: Path, output: Path, param1: float) -> None:
    tree = ET.parse(template)
    root = tree.getroot()
    nodes = root.findall(".//noise_level")
    if not nodes:
        raise RuntimeError(f"No <noise_level> node found in {template}")
    for node in nodes:
        node.set("param1", f"{param1:.12g}")
    output.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output, encoding="utf-8", xml_declaration=True)


def run_command(command: list[str], cwd: Path) -> None:
    print("[run]", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def format_command(template: str, run_dir: Path, sim_root: Path, sim_hittree: Path) -> str:
    return template.format(
        run_dir=str(run_dir),
        sim_root=str(sim_root),
        sim_hittree=str(sim_hittree),
    )


def prepare_or_run_simulation(cfg: ScanConfig, param1: float) -> Path:
    run_dir = cfg.outdir / "runs" / run_label(param1)
    run_dir.mkdir(parents=True, exist_ok=True)

    detector_parameters = run_dir / "detector_parameters.xml"
    write_detector_parameters(cfg.detector_parameters_template, detector_parameters, param1)

    sim_root = run_dir / cfg.sim_output_name
    sim_hittree = run_dir / cfg.sim_hittree_name

    if cfg.execute and not sim_root.exists():
        run_command(
            [
                cfg.ruby_command,
                str(cfg.simulation_runner),
                str(cfg.num_events),
                f"{cfg.energy_kev:g}",
                str(detector_parameters),
                str(cfg.random_seed),
                str(sim_root),
            ],
            cfg.simulation_runner.parent,
        )

    if cfg.execute and cfg.postprocess_command:
        command = format_command(cfg.postprocess_command, run_dir, sim_root, sim_hittree)
        print("[postprocess]", command)
        subprocess.run(command, cwd=run_dir, shell=True, check=True)

    if not sim_hittree.exists():
        print(f"[pending] {sim_hittree} does not exist yet.")
    return sim_hittree


def read_energies(root_path: Path, cfg: ScanConfig) -> np.ndarray:
    with uproot.open(root_path) as f:
        tree = f[cfg.tree_name]
        branches = [cfg.energy_branch, cfg.event_branch]
        arrays = tree.arrays(branches, library="np")

    energy = np.asarray(arrays[cfg.energy_branch], dtype=float)
    event = np.asarray(arrays[cfg.event_branch])

    if cfg.spectrum_mode == "hit":
        return energy
    if cfg.spectrum_mode != "event-sum":
        raise ValueError("spectrum.mode must be event-sum or hit.")

    order = np.argsort(event)
    event_sorted = event[order]
    energy_sorted = energy[order]
    unique_event, start, counts = np.unique(event_sorted,
                                            return_index=True,
                                            return_counts=True)
    del unique_event
    summed = np.add.reduceat(energy_sorted, start)
    if cfg.one_hit_only:
        summed = summed[counts == 1]
    return summed


def histogram(energies: np.ndarray, cfg: ScanConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.arange(cfg.histogram_min,
                      cfg.histogram_max + cfg.bin_width,
                      cfg.bin_width)
    counts, edges = np.histogram(energies, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return counts.astype(float), edges, centers


def compare_spectra(data_counts: np.ndarray,
                    sim_counts: np.ndarray,
                    centers: np.ndarray,
                    cfg: ScanConfig) -> tuple[float, float, int]:
    norm_mask = centers >= cfg.normalize_min
    data_integral = float(np.sum(data_counts[norm_mask]))
    sim_integral = float(np.sum(sim_counts[norm_mask]))
    alpha = data_integral / sim_integral if sim_integral > 0.0 else float("nan")

    fit_mask = (centers >= cfg.compare_min) & (centers <= cfg.compare_max)
    fit_mask &= np.isfinite(data_counts) & np.isfinite(sim_counts)
    if not np.isfinite(alpha):
        return alpha, float("nan"), 0

    residual = data_counts[fit_mask] - alpha * sim_counts[fit_mask]
    variance = np.maximum(data_counts[fit_mask], 1.0)
    chi2 = float(np.sum(residual * residual / variance))
    ndf = int(max(1, np.count_nonzero(fit_mask) - 1))
    return alpha, chi2, ndf


def write_scan_csv(results: list[ScanResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["param1", "chi2", "ndf", "chi2_ndf", "alpha",
                         "n_data", "n_sim", "sim_hittree"])
        for result in results:
            writer.writerow([
                result.param1,
                result.chi2,
                result.ndf,
                result.chi2_ndf,
                result.alpha,
                result.n_data,
                result.n_sim,
                result.sim_hittree,
            ])


def plot_scan(results: list[ScanResult], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = [result for result in results if np.isfinite(result.chi2_ndf)]
    if not valid:
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot([r.param1 for r in valid], [r.chi2_ndf for r in valid],
            marker="o", ls="-")
    ax.set_xlabel("noise_level param1")
    ax.set_ylabel("reduced chi-square")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_best_comparison(data_counts: np.ndarray,
                         sim_counts: np.ndarray,
                         edges: np.ndarray,
                         best: ScanResult,
                         path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.stairs(data_counts, edges, label="Data", color="black")
    ax.stairs(best.alpha * sim_counts, edges,
              label=f"Simulation x {best.alpha:.4g}, param1={best.param1:.5g}",
              color="tab:orange")
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    cfg = load_config(parse_args().config)
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    sim_hittrees = [
        (param1, prepare_or_run_simulation(cfg, param1))
        for param1 in cfg.param1_values
    ]
    sim_hittrees = [(param1, path) for param1, path in sim_hittrees if path.exists()]
    scan_csv = cfg.outdir / "scan_results.csv"
    if not sim_hittrees:
        write_scan_csv([], scan_csv)
        print(f"scan: {scan_csv}")
        print("Prepared detector_parameters.xml files, but no completed simulation hittrees were compared.")
        return

    data_energies = read_energies(cfg.data_hittree, cfg)
    data_counts, edges, centers = histogram(data_energies, cfg)
    results: list[ScanResult] = []
    best_sim_counts: np.ndarray | None = None
    best: ScanResult | None = None

    for param1, sim_hittree in sim_hittrees:
        sim_energies = read_energies(sim_hittree, cfg)
        sim_counts, _, _ = histogram(sim_energies, cfg)
        alpha, chi2, ndf = compare_spectra(data_counts, sim_counts, centers, cfg)
        result = ScanResult(param1=param1,
                            sim_hittree=sim_hittree,
                            alpha=alpha,
                            chi2=chi2,
                            ndf=ndf,
                            n_data=int(data_energies.size),
                            n_sim=int(sim_energies.size))
        results.append(result)
        print(f"param1={param1:.6g} chi2/ndf={result.chi2_ndf:.6g} alpha={alpha:.6g}")

        if np.isfinite(result.chi2_ndf) and (best is None or result.chi2_ndf < best.chi2_ndf):
            best = result
            best_sim_counts = sim_counts

    write_scan_csv(results, scan_csv)
    plot_scan(results, cfg.outdir / "chi2_scan.png")

    if best is not None and best_sim_counts is not None:
        plot_best_comparison(data_counts, best_sim_counts, edges,
                             best, cfg.outdir / "best_comparison.png")
        with (cfg.outdir / "best_summary.txt").open("w") as f:
            f.write(f"best_param1: {best.param1:.12g}\n")
            f.write(f"chi2: {best.chi2:.12g}\n")
            f.write(f"ndf: {best.ndf}\n")
            f.write(f"chi2_ndf: {best.chi2_ndf:.12g}\n")
            f.write(f"alpha: {best.alpha:.12g}\n")
            f.write(f"sim_hittree: {best.sim_hittree}\n")

    print(f"scan: {scan_csv}")
    if best is None:
        print("No completed simulation hittrees were compared. "
              "Set simulation.execute: true or place hittrees in the run directories.")
    else:
        print(f"best param1: {best.param1:.6g} (chi2/ndf={best.chi2_ndf:.6g})")


if __name__ == "__main__":
    main()
