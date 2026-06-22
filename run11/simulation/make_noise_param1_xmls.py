#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import csv
import xml.etree.ElementTree as ET


def xml_stylesheet_line(path: Path) -> str:
    for line in path.read_text().splitlines():
        if line.strip().startswith("<?xml-stylesheet"):
            return line
    return ""


def label(param1: float) -> str:
    return f"param1_{param1:.3f}".replace(".", "p")


def write_xml(template_xml: Path, output_xml: Path, param1: float) -> None:
    tree = ET.parse(template_xml)
    root = tree.getroot()
    for node in root.findall(".//noise_level"):
        node.set("param1", f"{param1:.12g}")
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(tree, space="  ")

    lines = ['<?xml version="1.0" ?>']
    stylesheet = xml_stylesheet_line(template_xml)
    if stylesheet:
        lines.append(stylesheet)
    lines.append(ET.tostring(root, encoding="unicode"))
    output_xml.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    simulation_dir = Path(__file__).resolve().parent
    run11_dir = simulation_dir.parent

    template_xml = run11_dir / "database/detector_parameters/detector_parameters.xml"
    output_root = simulation_dir / "products/noise_param1_scan"

    #param1_start = 0.020
    #param1_stop  = 0.080
    param1_start = 0.022
    param1_stop  = 0.080
    param1_step  = 0.002
    epsilon      = 1e-5

    param1_values = np.arange(param1_start, param1_stop + epsilon, param1_step)

    rows = []
    for param1 in param1_values:
        run_dir = output_root / label(param1)
        output_xml = run_dir / "detector_parameters.xml"
        write_xml(template_xml, output_xml, param1)
        rows.append([param1, output_xml])
        print(f"{param1:.3f} -> {output_xml}")

    with (output_root / "param1_files.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["param1", "detector_parameters_xml"])
        writer.writerows(rows)
