#!/usr/bin/env python3
from pathlib import Path
import csv

import numpy as np
import uproot


def load_energies(root_path: Path, tree_name: str, energy_branch: str,
                  one_hit_only: bool, section: int | None) -> np.ndarray:
    with uproot.open(root_path) as f:
        tree = f[tree_name]
        branches = [energy_branch]
        if one_hit_only:
            branches.append("num_hits")
        if section is not None:
            branches.append("section")
        arrays = tree.arrays(branches, library="np")

    energy = np.asarray(arrays[energy_branch], dtype=float)
    mask = np.isfinite(energy)
    if one_hit_only:
        mask &= arrays["num_hits"] == 1
    if section is not None:
        mask &= arrays["section"] == section
    return energy[mask]


def hist(energies: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.histogram(energies, bins=edges)[0].astype(float)


def fit_scale(data_energy: np.ndarray, sim_counts: np.ndarray, edges: np.ndarray,
              scale_values: np.ndarray, norm_min: float, fit_min: float,
              fit_max: float) -> tuple[float, float, float, int]:
    centers   = 0.5 * (edges[:-1] + edges[1:])
    norm_mask = centers >= norm_min
    fit_mask  = (centers >= fit_min) & (centers <= fit_max)

    best = (float("nan"), float("nan"), float("nan"), 0)
    for scale in scale_values:
        data_counts = hist(data_energy * scale, edges)
        data_sum    = data_counts[norm_mask].sum()
        sim_sum     = sim_counts[norm_mask].sum()
        alpha       = data_sum / sim_sum if sim_sum > 0 else float("nan")
        if not np.isfinite(alpha):
            continue

        residual = data_counts[fit_mask] - alpha * sim_counts[fit_mask]
        variance = np.maximum(data_counts[fit_mask], 1.0)
        chi2     = float(np.sum(residual * residual / variance))
        ndf      = max(1, int(np.count_nonzero(fit_mask) - 2))
        chi2_ndf = chi2 / ndf
        if not np.isfinite(best[2]) or chi2_ndf < best[2]:
            best = (scale, alpha, chi2_ndf, ndf)
    return best


def plot_best(data_energy: np.ndarray, sim_counts: np.ndarray, edges: np.ndarray,
              scale: float, alpha: float, title: str, output_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data_counts = hist(data_energy * scale, edges)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.errorbar(0.5 * (edges[:-1] + edges[1:]), data_counts,
                yerr=np.sqrt(np.maximum(data_counts, 1.0)),
                fmt=".", ms=3, color="black", label="Data")
    ax.stairs(alpha * sim_counts, edges, color="tab:orange", label="Simulation")
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_png, dpi=160)
    plt.close(fig)


def param1_from_label(path: Path) -> float:
    return float(path.name.replace("param1_", "").replace("p", "."))


if __name__ == "__main__":
    simulation_dir = Path(__file__).resolve().parent
    run11_dir = simulation_dir.parent

    scan_root  = simulation_dir / "products/noise_param1_scan"
    data_root  = run11_dir / "products/hittree_merge_Na22.root"
    output_csv = scan_root / "fit_results.csv"
    output_png = scan_root / "best_fit.png"

    tree_name     = "hittree"
    energy_branch = "energy"
    one_hit_only  = True
    section       = None

    energy_min = 0.0
    energy_max = 1600.0
    bin_width  = 20.0

    scale_min  = 0.90
    scale_max  = 1.10
    scale_step = 0.001

    norm_min = 350.0
    fit_min  = 350.0
    fit_max  = 1300.0

    edges  = np.arange(energy_min, energy_max + bin_width, bin_width)
    scales = np.arange(scale_min, scale_max + 0.5 * scale_step, scale_step)
    data_energy = load_energies(data_root, tree_name, energy_branch, one_hit_only, section)

    rows = []
    best_row = None
    best_sim_counts = None
    sim_root_name = "simulation_Na22_merged.root"

    for run_dir in sorted(scan_root.glob("param1_*")):
        sim_root = run_dir / sim_root_name
        if not sim_root.exists():
            continue
        sim_energy = load_energies(sim_root, tree_name, energy_branch, one_hit_only, section)
        sim_counts = hist(sim_energy, edges)
        scale, alpha, chi2_ndf, ndf = fit_scale(
            data_energy, sim_counts, edges, scales, norm_min, fit_min, fit_max
        )
        param1 = param1_from_label(run_dir)
        rows.append([param1, scale, alpha, chi2_ndf, ndf, data_energy.size, sim_energy.size, sim_root])
        print(f"param1={param1:.3f} scale={scale:.4f} alpha={alpha:.4g} chi2/ndf={chi2_ndf:.4g}")
        if best_row is None or chi2_ndf < best_row[3]:
            best_row = rows[-1]
            best_sim_counts = sim_counts

    with output_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["param1", "scale", "alpha", "chi2_ndf", "ndf",
                         "n_data", "n_sim", "sim_root"])
        writer.writerows(rows)

    if best_row is not None:
        title = f"param1={best_row[0]:.3f}, scale={best_row[1]:.4f}, chi2/ndf={best_row[3]:.3g}"
        plot_best(data_energy, best_sim_counts, edges, best_row[1], best_row[2], title, output_png)
        print(f"best: {title}")
    print(f"saved: {output_csv}")
