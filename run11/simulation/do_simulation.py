#!/usr/bin/env python3
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import uproot


def run_simulation_energy(energy_str, num_photons, detector_parameter_file, output):
    output = Path(output)

    cmd_list = [
        "./run_simulation.rb",
        str(num_photons),
        energy_str,
        detector_parameter_file,
        str(output),
    ]

    print(f"Start: {energy_str} keV, {num_photons:.0f} photons, {detector_parameter_file}")

    result = subprocess.run(
        cmd_list,
        check=True,
    )

    print(f"Done : {energy_str} keV")

    return energy_str, output, result.returncode


def remap_eventid(arrays, next_eventid):
    arrays = dict(arrays)
    if "eventid" not in arrays or len(arrays["eventid"]) == 0:
        return arrays, next_eventid

    _, inverse = np.unique(arrays["eventid"], return_inverse=True)
    arrays["eventid"] = (inverse + next_eventid).astype(arrays["eventid"].dtype, copy=False)
    return arrays, next_eventid + int(inverse.max()) + 1


def merge_hittrees(input_files, output_file, tree_name="hittree"):
    merged_parts = []
    next_eventid = 0

    for input_file in input_files:
        with uproot.open(input_file) as f:
            arrays = f[tree_name].arrays(library="np")
        arrays, next_eventid = remap_eventid(arrays, next_eventid)
        merged_parts.append(arrays)

    branch_names = list(merged_parts[0].keys())
    merged = {
        name: np.concatenate([part[name] for part in merged_parts])
        for name in branch_names
    }
    branch_types = {name: values.dtype for name, values in merged.items()}

    with uproot.recreate(output_file) as f:
        f.mktree(tree_name, branch_types)
        f[tree_name].extend(merged)

    print(f"Merged: {output_file} ({len(next(iter(merged.values())))} hits)")


def build_jobs(detector_parameter_files):
    jobs = []
    merge_inputs = {}
    output_to_merge = {}

    for detector_parameter_file in detector_parameter_files:
        for source_name, energy_dict in simConfigDict.items():
            outputs = []
            for energy_str, num_photons in energy_dict.items():
                fileName = f"simulation_{source_name}_{int(float(energy_str))}keV.root"
                output = Path(detector_parameter_file).parent / fileName
                outputs.append(output)
                jobs.append((source_name, energy_str, int(num_photons), detector_parameter_file, output))

            merged_output = Path(detector_parameter_file).parent / f"simulation_{source_name}_merged.root"
            merge_inputs[merged_output] = outputs
            for output in outputs:
                output_to_merge[output] = merged_output

    return jobs, merge_inputs, output_to_merge


def run_all_param1(detector_parameter_files):
    jobs, merge_inputs, output_to_merge = build_jobs(detector_parameter_files)
    finished_outputs = {merged_output: set() for merged_output in merge_inputs}
    merged_outputs = set()

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {}
        for source_name, energy_str, num_photons, detector_parameter_file, output in jobs:
            print(output)
            future = executor.submit(run_simulation_energy, energy_str, num_photons, detector_parameter_file, output)
            futures[future] = output

        for future in as_completed(futures):
            energy_str, output, returncode = future.result()
            print(f"Finished: {energy_str} keV, returncode={returncode}")

            merged_output = output_to_merge[output]
            finished_outputs[merged_output].add(output)
            if (
                merged_output not in merged_outputs
                and len(finished_outputs[merged_output]) == len(merge_inputs[merged_output])
            ):
                merge_hittrees(merge_inputs[merged_output], merged_output, tree_name=tree_name)
                merged_outputs.add(merged_output)


if __name__ == "__main__":

    num_workers = 10
    tree_name = "hittree"

    num_photons = 1e9

    simConfigDict = {
        "Na22": {"511.0": 2*num_photons,  
                 "1275.4": 1*num_photons}
    }

    metaFile = "products/noise_param1_scan/param1_files.csv"
    with open(metaFile, "r") as f:
        lines = f.readlines()
        num_rows = len(lines)
        detector_parameter_files = []

        for iline in range(1, num_rows):
            line_list = lines[iline].replace("\n","").split(",")
            _, detector_parameter_file = line_list
            print(detector_parameter_file)
            detector_parameter_files.append(detector_parameter_file)

    run_all_param1(detector_parameter_files)
