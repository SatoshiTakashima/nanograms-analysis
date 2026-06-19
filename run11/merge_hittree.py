#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import uproot


def read_time_ids(csv_path: Path, time_column: str) -> list[str]:
    table = pd.read_csv(csv_path)
    return table[time_column].astype(str).tolist()


def input_path(products_root: Path, time_id: str, hittree_name: str) -> Path:
    return products_root / time_id / hittree_name


def remap_eventid(eventid: np.ndarray, next_eventid: int) -> tuple[np.ndarray, int]:
    """Map event IDs in one input file to new IDs starting at next_eventid."""
    unique_eventids, inverse = np.unique(eventid, return_inverse=True)
    new_eventid = inverse.astype(eventid.dtype, copy=False) + np.asarray(next_eventid, dtype=eventid.dtype)
    return new_eventid, next_eventid + len(unique_eventids)


def write_ttree(output_path: Path, tree_name: str, arrays: dict[str, np.ndarray]) -> None:
    branch_types = {name: values.dtype for name, values in arrays.items()}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as output_file:
        tree = output_file.mktree(tree_name, branch_types)
        tree.extend(arrays)


def merge_hittrees(
    *,
    time_ids: list[str],
    products_root: Path,
    hittree_name: str,
    tree_name: str,
    eventid_branch: str,
    output_path: Path,
) -> tuple[int, int]:
    first_path = input_path(products_root, time_ids[0], hittree_name)
    with uproot.open(first_path) as first_file:
        first_tree = first_file[tree_name]
        branch_names = first_tree.keys()
        empty_arrays = {
            name: first_tree[name].array(library="np", entry_stop=0)
            for name in branch_names
        }

    pieces = {name: [] for name in branch_names}
    total_entries = 0
    next_eventid  = 0

    for time_id in time_ids:
        path = input_path(products_root, time_id, hittree_name)
        if not path.exists():
            raise FileNotFoundError(path)

        with uproot.open(path) as root_file:
            tree = root_file[tree_name]
            arrays = tree.arrays(branch_names, library="np")

        remapped_eventid, next_eventid = remap_eventid(arrays[eventid_branch], next_eventid)
        for name in branch_names:
            pieces[name].append(remapped_eventid if name == eventid_branch else arrays[name])
        total_entries += len(arrays[eventid_branch])

    merged = {
        name: np.concatenate(values) if values else empty_arrays[name]
        for name, values in pieces.items()
    }
    write_ttree(output_path, tree_name, merged)

    return total_entries, next_eventid


if __name__ == "__main__":
    run11_dir     = Path(__file__).resolve().parent
    csv_path      = run11_dir / "metadata/data_group/data_group_Na22.csv"
    products_root = run11_dir / "products"
    output_path   = run11_dir / "products/hittree_merge_Na22.root"

    time_ids = read_time_ids(csv_path, "time")
    entries, events = merge_hittrees(
        time_ids=time_ids,
        products_root=products_root,
        hittree_name="hittree.root",
        tree_name="hittree",
        eventid_branch="eventid",
        output_path=output_path,
    )

    print(f"{entries} entries, {events} events -> {output_path} (wrote)")
