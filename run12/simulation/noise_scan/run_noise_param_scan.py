#!/usr/bin/env python3
"""Scan detector noise_level param1 and compare simulation hittrees to data.

For each configured param1 value this script:
  1. copies detector_parameters.xml and edits <noise_level param1="...">,
  2. writes a small Ruby simulation runner using that modified XML,
  3. optionally executes the simulation and an optional postprocess command,
  4. reads the produced simulation hittree,
  5. compares the 1-hit energy histogram to the experimental hittree,
  6. saves reduced chi-square results and best-fit plots.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import uproot
import yaml


@dataclass(frozen=True)
class ScanConfig:
    outdir: Path
    execute: bool
    num_events: int
    random_seed: int
    ruby_command: str
    detector_configuration: Path
    detector_parameters_template: Path
    gdml: Path
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


def load_config(path: str | Path) -> ScanConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as f:
        config = yaml.safe_load(f) or {}

    scan_node = nested_get(config, "scan", {}) or {}
    if "values" in scan_node:
        param1_values = [float(value) for value in scan_node["values"]]
    else:
        param1_values = np.linspace(
            float(scan_node.get("min", 0.02)),
            float(scan_node.get("max", 0.08)),
            int(scan_node.get("steps", 7)),
        ).tolist()

    return ScanConfig(
        outdir=resolve_path(config_path, nested_get(config, "output.outdir", "products/noise_scan")),
        execute=bool(nested_get(config, "simulation.execute", False)),
        num_events=int(nested_get(config, "simulation.num_events", 10000)),
        random_seed=int(nested_get(config, "simulation.random_seed", 0)),
        ruby_command=str(nested_get(config, "simulation.ruby_command", "ruby")),
        detector_configuration=resolve_path(
            config_path,
            required_value(config, "simulation.detector_configuration"),
        ),
        detector_parameters_template=resolve_path(
            config_path,
            required_value(config, "simulation.detector_parameters_template"),
        ),
        gdml=resolve_path(config_path, required_value(config, "simulation.gdml")),
        param1_values=param1_values,
        sim_output_name=str(nested_get(config, "simulation.output_root", "simulation.root")),
        sim_hittree_name=str(nested_get(config, "simulation.hittree_root", "simulation.root")),
        postprocess_command=nested_get(config, "simulation.postprocess_command", None),
        data_hittree=resolve_path(config_path, required_value(config, "data.hittree")),
        tree_name=str(nested_get(config, "data.tree", "hittree")),
        energy_branch=str(nested_get(config, "branches.energy", "energy")),
        event_branch=str(nested_get(config, "branches.event", "eventid")),
        histogram_min=float(nested_get(config, "histogram.min_kev", 0.0)),
        histogram_max=float(nested_get(config, "histogram.max_kev", 2000.0)),
        bin_width=float(nested_get(config, "histogram.bin_width_kev", 5.0)),
        compare_min=float(nested_get(config, "comparison.min_kev", 350.0)),
        compare_max=float(nested_get(config, "comparison.max_kev", 1300.0)),
        normalize_min=float(nested_get(config, "comparison.normalize_min_kev", 350.0)),
        spectrum_mode=str(nested_get(config, "spectrum.mode", "event-sum")),
        one_hit_only=bool(nested_get(config, "spectrum.one_hit_only", True)),
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


def write_ruby_runner(cfg: ScanConfig, run_dir: Path, detector_parameters: Path) -> Path:
    script = run_dir / "run_simulation.rb"
    sim_output = run_dir / cfg.sim_output_name
    content = f"""#! /usr/bin/env ruby

require 'comptonsoft'

num = {cfg.num_events}
random = {cfg.random_seed}
energy = 1332.5

sim = ComptonSoft::Simulation.new
sim.output = {str(sim_output)!r}
sim.random_seed = random
sim.verbose = 0
sim.set_database(detector_configuration: {str(cfg.detector_configuration)!r},
                 detector_parameters: {str(detector_parameters)!r})
sim.set_gdml {str(cfg.gdml)!r}

sim.set_physics(hadron_hp: false, cut_value: 0.001, radioactive_decay: true)
sim.enable_timing_process
sim.set_primary_generator :PointSourcePrimaryGen, {{
  particle: "gamma",
  spectral_distribution: "gaussian",
  energy_mean: energy,
  energy_sigma: energy*0.001,
  position: vec(159.4, 0.0, 37.5),
  direction: vec(-1.0, 0.0, 0.0),
  theta_min: 0.0,
  theta_max: 90.0*Math::PI/180.0,
}}

sim.run(num)
"""
    script.write_text(content)
    script.chmod(0o755)
    return script


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
    ruby_script = write_ruby_runner(cfg, run_dir, detector_parameters)

    sim_root = run_dir / cfg.sim_output_name
    sim_hittree = run_dir / cfg.sim_hittree_name

    if cfg.execute and not sim_root.exists():
        run_command([cfg.ruby_command, str(ruby_script)], run_dir)

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

    data_energies = read_energies(cfg.data_hittree, cfg)
    data_counts, edges, centers = histogram(data_energies, cfg)

    results: list[ScanResult] = []
    best_sim_counts: np.ndarray | None = None
    best: ScanResult | None = None

    for param1 in cfg.param1_values:
        sim_hittree = prepare_or_run_simulation(cfg, param1)
        if not sim_hittree.exists():
            continue

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

    scan_csv = cfg.outdir / "scan_results.csv"
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
