#!/usr/bin/env python3
"""Split a tpctree ROOT file into chunks small enough for quicklook output."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import uproot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split tpctree ROOT files.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--max-entries", type=int, default=4000)
    parser.add_argument("--config", type=Path, default=None)
    return parser.parse_args()


def chunk_stops(trigger_ids, max_entries: int) -> list[tuple[int, int]]:
    n_entries = len(trigger_ids)
    chunks: list[tuple[int, int]] = []
    start = 0
    while start < n_entries:
        stop = min(start + max_entries, n_entries)
        if stop < n_entries:
            while stop > start + 1 and trigger_ids[stop - 1] == trigger_ids[stop]:
                stop -= 1
        chunks.append((start, stop))
        start = stop
    return chunks


def root_clone_tree(input_path: Path, output_path: Path, start: int, stop: int) -> None:
    macro_path = output_path.with_suffix(".C")
    macro_path.write_text(
        "\n".join(
            [
                "{",
                f"  TFile input({json.dumps(str(input_path))}, \"READ\");",
                '  TTree* tree = nullptr;',
                '  input.GetObject("tpctree", tree);',
                '  if (!tree) { throw std::runtime_error("missing tpctree"); }',
                f"  TFile output({json.dumps(str(output_path))}, \"RECREATE\");",
                "  TTree* chunk = tree->CloneTree(0);",
                f"  for (Long64_t entry = {start}; entry < {stop}; ++entry) {{",
                "    tree->GetEntry(entry);",
                "    chunk->Fill();",
                "  }",
                '  chunk->Write("tpctree");',
                "  output.Close();",
                "  input.Close();",
                "}",
                "",
            ]
        )
    )
    try:
        subprocess.run(
            ["root", "-l", "-b", "-q", str(macro_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(exc.stdout) from exc
    finally:
        macro_path.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with uproot.open(args.input) as root_file:
        tree = root_file["tpctree"]
        n_entries = tree.num_entries
        if n_entries <= args.max_entries:
            print(f"000\t{args.input}")
            return

        trigger_ids = tree["triggerid"].array(library="np")
        chunks = chunk_stops(trigger_ids, args.max_entries)
        for index, (start, stop) in enumerate(chunks):
            part_dir = args.output_dir / f"part{index:03d}"
            part_dir.mkdir(parents=True, exist_ok=True)
            if args.config and args.config.exists():
                shutil.copy2(args.config, part_dir / args.config.name)
            output = part_dir / args.input.name
            root_clone_tree(args.input, output, start, stop)
            print(f"{index:03d}\t{output}")


if __name__ == "__main__":
    main()
