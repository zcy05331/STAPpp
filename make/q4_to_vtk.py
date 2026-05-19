#!/usr/bin/env python3
"""Convert one STAPpp Q4 .dat/.out pair to legacy VTK for ParaView."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple


def parse_dat(path: Path):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    i = 0
    i += 1  # title line may contain spaces
    numnp = int(lines[i][0]); numeg = int(lines[i][1]); nlcase = int(lines[i][2]); modex = int(lines[i][3]); i += 1
    nodes = []
    for _ in range(numnp):
        parts = lines[i]
        node = int(parts[0]); x = float(parts[4]); y = float(parts[5]); z = float(parts[6]); i += 1
        nodes.append((node, x, y, z))
    for _ in range(nlcase):
        i += 1  # load case number
        nloads = int(lines[i][0]); i += 1 + nloads
    elements = []
    etype = int(lines[i][0]); nume = int(lines[i][1]); nummat = int(lines[i][2]); i += 1
    if etype != 2:
        raise SystemExit(f"expected Q4 element type 2, got {etype}")
    i += nummat
    for _ in range(nume):
        parts = lines[i]
        ele = int(parts[0]); n1 = int(parts[1]); n2 = int(parts[2]); n3 = int(parts[3]); n4 = int(parts[4]); i += 1
        elements.append((ele, n1, n2, n3, n4))
    return nodes, elements


def parse_out(path: Path):
    disps: Dict[int, Tuple[float, float, float]] = {}
    stresses: Dict[int, List[Tuple[float, float, float]]] = {}
    in_disp = False
    for line in path.read_text().splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_disp = True
            continue
        if in_disp and "S T R E S S" in line:
            in_disp = False
        parts = line.split()
        if in_disp and len(parts) == 4 and parts[0].isdigit():
            disps[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
        if len(parts) == 5 and parts[0].isdigit() and parts[1].isdigit():
            ele = int(parts[0])
            stresses.setdefault(ele, []).append((float(parts[2]), float(parts[3]), float(parts[4])))
    return disps, stresses


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dat", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("vtk", type=Path)
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()

    nodes, elements = parse_dat(args.dat)
    disps, stresses = parse_out(args.out)
    args.vtk.parent.mkdir(parents=True, exist_ok=True)

    with args.vtk.open("w") as f:
        f.write("# vtk DataFile Version 3.0\nSTAPpp Q4 result\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(nodes)} double\n")
        for node, x, y, z in nodes:
            ux, uy, uz = disps.get(node, (0.0, 0.0, 0.0))
            f.write(f"{x + args.scale * ux:.12e} {y + args.scale * uy:.12e} {z + args.scale * uz:.12e}\n")
        f.write(f"CELLS {len(elements)} {len(elements) * 5}\n")
        for _ele, n1, n2, n3, n4 in elements:
            f.write(f"4 {n1-1} {n2-1} {n3-1} {n4-1}\n")
        f.write(f"CELL_TYPES {len(elements)}\n")
        for _ in elements:
            f.write("9\n")
        f.write(f"POINT_DATA {len(nodes)}\nVECTORS displacement double\n")
        for node, *_ in nodes:
            ux, uy, uz = disps.get(node, (0.0, 0.0, 0.0))
            f.write(f"{ux:.12e} {uy:.12e} {uz:.12e}\n")
        f.write(f"CELL_DATA {len(elements)}\nSCALARS sigma_x double 1\nLOOKUP_TABLE default\n")
        for ele, *_ in elements:
            vals = stresses.get(ele, [(0.0, 0.0, 0.0)])
            f.write(f"{sum(v[0] for v in vals)/len(vals):.12e}\n")
        f.write("SCALARS von_mises_plane_stress double 1\nLOOKUP_TABLE default\n")
        for ele, *_ in elements:
            vals = stresses.get(ele, [(0.0, 0.0, 0.0)])
            sx = sum(v[0] for v in vals)/len(vals)
            sy = sum(v[1] for v in vals)/len(vals)
            txy = sum(v[2] for v in vals)/len(vals)
            vm = (sx*sx - sx*sy + sy*sy + 3.0*txy*txy) ** 0.5
            f.write(f"{vm:.12e}\n")
    print(f"wrote {args.vtk}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
