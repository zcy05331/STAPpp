#!/usr/bin/env python3
"""Generate and validate STAPpp Q4 course-project cases using only stdlib."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
E = 1000.0
NU = 0.25
THICKNESS = 1.0
SIGMA_X = 1.0

Node = Tuple[int, int, int, int, float, float, float]
Load = Tuple[int, int, float]
Element = Tuple[int, int, int, int, int, int]


def rel_abs_error(value: float, expected: float) -> float:
    denom = max(1.0, abs(expected))
    return abs(value - expected) / denom


def write_case(case_dir: Path, name: str, title: str, nodes: List[Node], loads: List[Load], elements: List[Element]) -> Path:
    case_dir.mkdir(parents=True, exist_ok=True)
    dat = case_dir / f"{name}.dat"
    lines: List[str] = [title]
    lines.append(f"{len(nodes)} 1 1 1")
    for node, bx, by, bz, x, y, z in nodes:
        lines.append(f"{node:5d} {bx:2d} {by:2d} {bz:2d} {x:18.10e} {y:18.10e} {z:18.10e}")
    lines.append("1")
    lines.append(str(len(loads)))
    for node, dof, load in loads:
        lines.append(f"{node:5d} {dof:2d} {load:18.10e}")
    lines.append(f"2 {len(elements)} 1")
    lines.append(f"1 {E:.10e} {NU:.10e} {THICKNESS:.10e} 0")
    for ele, n1, n2, n3, n4, mset in elements:
        lines.append(f"{ele:5d} {n1:5d} {n2:5d} {n3:5d} {n4:5d} {mset:5d}")
    dat.write_text("\n".join(lines) + "\n")
    return dat


def generate_rect_mesh(nx: int, ny: int, *, L: float, H: float, patch_bc: bool, vertical_load: bool = False) -> Tuple[List[Node], List[Element], Dict[Tuple[int, int], int]]:
    nodes: List[Node] = []
    node_id: Dict[Tuple[int, int], int] = {}
    n = 1
    for j in range(ny + 1):
        y = H * j / ny
        for i in range(nx + 1):
            x = L * i / nx
            if patch_bc:
                bx = 1 if i == 0 else 0
                by = 1 if (i == 0 and j == 0) else 0
            else:
                bx = 1 if i == 0 else 0
                by = 1 if i == 0 else 0
            bz = 1
            node_id[(i, j)] = n
            nodes.append((n, bx, by, bz, x, y, 0.0))
            n += 1

    elements: List[Element] = []
    e = 1
    for j in range(ny):
        for i in range(nx):
            n1 = node_id[(i, j)]
            n2 = node_id[(i + 1, j)]
            n3 = node_id[(i + 1, j + 1)]
            n4 = node_id[(i, j + 1)]
            elements.append((e, n1, n2, n3, n4, 1))
            e += 1
    return nodes, elements, node_id


def edge_loads_x(node_id: Dict[Tuple[int, int], int], nx: int, ny: int, H: float, total: float) -> List[Load]:
    loads = {node_id[(nx, j)]: 0.0 for j in range(ny + 1)}
    traction = total / H
    for j in range(ny):
        length = H / ny
        f = traction * length * THICKNESS / 2.0
        loads[node_id[(nx, j)]] += f
        loads[node_id[(nx, j + 1)]] += f
    return [(node, 1, load) for node, load in sorted(loads.items())]


def edge_loads_y(node_id: Dict[Tuple[int, int], int], nx: int, ny: int, H: float, total: float) -> List[Load]:
    loads = {node_id[(nx, j)]: 0.0 for j in range(ny + 1)}
    traction = total / H
    for j in range(ny):
        length = H / ny
        f = traction * length * THICKNESS / 2.0
        loads[node_id[(nx, j)]] += f
        loads[node_id[(nx, j + 1)]] += f
    return [(node, 2, load) for node, load in sorted(loads.items())]


def run_case(exe: Path, dat: Path) -> Path:
    base = dat.with_suffix("")
    out = dat.with_suffix(".out")
    if out.exists():
        out.unlink()
    result = subprocess.run([str(exe), str(base)], cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"case {dat} failed with code {result.returncode}\n{result.stdout}")
    if not out.exists():
        raise RuntimeError(f"case {dat} did not produce {out}")
    return out


def parse_displacements(out: Path) -> Dict[int, Tuple[float, float, float]]:
    disp: Dict[int, Tuple[float, float, float]] = {}
    in_section = False
    for line in out.read_text().splitlines():
        if "D I S P L A C E M E N T S" in line:
            in_section = True
            continue
        if in_section and "S T R E S S" in line:
            break
        parts = line.split()
        if in_section and len(parts) == 4 and parts[0].isdigit():
            disp[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
    return disp


def parse_q4_stresses(out: Path) -> List[Tuple[int, int, float, float, float]]:
    rows: List[Tuple[int, int, float, float, float]] = []
    for line in out.read_text().splitlines():
        parts = line.split()
        if len(parts) == 5 and parts[0].isdigit() and parts[1].isdigit():
            rows.append((int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
    return rows


def check_patch(exe: Path, case_dir: Path, name: str, nx: int, ny: int, stress_tol: float, disp_tol: float) -> Dict[str, float]:
    nodes, elements, node_id = generate_rect_mesh(nx, ny, L=1.0, H=1.0, patch_bc=True)
    loads = edge_loads_x(node_id, nx, ny, H=1.0, total=SIGMA_X)
    dat = write_case(case_dir, name, f"Q4 {name} constant stress patch", nodes, loads, elements)
    out = run_case(exe, dat)
    stresses = parse_q4_stresses(out)
    if len(stresses) != len(elements) * 4:
        raise AssertionError(f"{name}: expected {len(elements) * 4} stress rows, got {len(stresses)}")
    max_stress_err = 0.0
    for _, _, sx, sy, txy in stresses:
        max_stress_err = max(max_stress_err, rel_abs_error(sx, SIGMA_X), abs(sy), abs(txy))
    if max_stress_err > stress_tol:
        raise AssertionError(f"{name}: stress error {max_stress_err:.3e} > {stress_tol:.3e}")

    disp = parse_displacements(out)
    max_disp_err = 0.0
    for node, _, _, _, x, y, _ in nodes:
        ux_expected = SIGMA_X / E * x
        uy_expected = -NU * SIGMA_X / E * y
        ux, uy, _ = disp[node]
        max_disp_err = max(max_disp_err, rel_abs_error(ux, ux_expected), rel_abs_error(uy, uy_expected))
    if max_disp_err > disp_tol:
        raise AssertionError(f"{name}: displacement error {max_disp_err:.3e} > {disp_tol:.3e}")

    return {"max_stress_error": max_stress_err, "max_displacement_error": max_disp_err}


def check_skew_patch(exe: Path) -> Dict[str, float]:
    case_dir = ROOT / "data/q4_patch_skew"
    nodes: List[Node] = [
        (1, 1, 1, 1, 0.0, 0.0, 0.0),
        (2, 0, 0, 1, 2.0, 0.0, 0.0),
        (3, 0, 0, 1, 2.5, 1.0, 0.0),
        (4, 1, 0, 1, 0.0, 1.0, 0.0),
    ]
    elements: List[Element] = [(1, 1, 2, 3, 4, 1)]
    loads: List[Load] = [
        (2, 1, 0.5 * SIGMA_X * THICKNESS),
        (3, 1, 0.5 * SIGMA_X * THICKNESS),
    ]
    dat = write_case(case_dir, "q4_patch_skew", "Q4 skew quadrilateral constant stress patch", nodes, loads, elements)
    out = run_case(exe, dat)

    stresses = parse_q4_stresses(out)
    if len(stresses) != 4:
        raise AssertionError(f"q4_patch_skew: expected 4 stress rows, got {len(stresses)}")

    max_stress_err = 0.0
    for _, _, sx, sy, txy in stresses:
        max_stress_err = max(max_stress_err, rel_abs_error(sx, SIGMA_X), abs(sy), abs(txy))
    if max_stress_err > 1e-8:
        raise AssertionError(f"q4_patch_skew: stress error {max_stress_err:.3e} > 1e-8")

    disp = parse_displacements(out)
    max_disp_err = 0.0
    for node, _, _, _, x, y, _ in nodes:
        ux_expected = SIGMA_X / E * x
        uy_expected = -NU * SIGMA_X / E * y
        ux, uy, _ = disp[node]
        max_disp_err = max(max_disp_err, rel_abs_error(ux, ux_expected), rel_abs_error(uy, uy_expected))
    if max_disp_err > 1e-8:
        raise AssertionError(f"q4_patch_skew: displacement error {max_disp_err:.3e} > 1e-8")

    return {"max_stress_error": max_stress_err, "max_displacement_error": max_disp_err}


def run_convergence(exe: Path) -> List[Tuple[int, int, float, float, float]]:
    case_dir = ROOT / "data/q4_convergence"
    L = 4.0
    H = 1.0
    meshes = [(2, 1), (4, 2), (8, 2), (16, 4), (32, 8)]
    values: List[Tuple[int, int, float]] = []
    for nx, ny in meshes:
        nodes, elements, node_id = generate_rect_mesh(nx, ny, L=L, H=H, patch_bc=False)
        loads = edge_loads_y(node_id, nx, ny, H=H, total=-1.0)
        name = f"q4_cantilever_{nx}x{ny}"
        dat = write_case(case_dir, name, f"Q4 cantilever convergence {nx}x{ny}", nodes, loads, elements)
        out = run_case(exe, dat)
        disp = parse_displacements(out)
        tip = node_id[(nx, ny)]
        values.append((nx, ny, disp[tip][1]))

    ref = values[-1][2]
    table: List[Tuple[int, int, float, float, float]] = []
    prev_err = None
    prev_h = None
    for nx, ny, uy in values[:-1]:
        err = abs(uy - ref) / max(abs(ref), 1.0)
        h = L / nx
        order = float("nan")
        if prev_err is not None and err > 0 and prev_err > 0:
            order = math.log(prev_err / err) / math.log(prev_h / h)
        table.append((nx, ny, uy, err, order))
        prev_err = err
        prev_h = h

    errors = [row[3] for row in table]
    if any(errors[i + 1] > errors[i] * 1.02 for i in range(len(errors) - 1)):
        raise AssertionError(f"convergence errors are not monotonic enough: {errors}")
    return table



def expect_failure(exe: Path, dat: Path, expected_text: str) -> None:
    base = dat.with_suffix("")
    result = subprocess.run([str(exe), str(base)], cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode == 0:
        raise AssertionError(f"negative case {dat} unexpectedly succeeded")
    if expected_text not in result.stdout:
        raise AssertionError(f"negative case {dat} did not include expected text {expected_text!r}\n{result.stdout}")


def run_invalid_cases(exe: Path) -> None:
    case_dir = ROOT / "data/q4_invalid"
    case_dir.mkdir(parents=True, exist_ok=True)

    nodes, elements, node_id = generate_rect_mesh(1, 1, L=1.0, H=1.0, patch_bc=True)
    loads = edge_loads_x(node_id, 1, 1, H=1.0, total=SIGMA_X)

    active_z_nodes = list(nodes)
    n, bx, by, _bz, x, y, z = active_z_nodes[1]
    active_z_nodes[1] = (n, bx, by, 0, x, y, z)
    dat = write_case(case_dir, "q4_invalid_active_z", "Q4 invalid active z DOF", active_z_nodes, loads, elements)
    expect_failure(exe, dat, "Q4 element requires the z DOF")

    bad_elements = [(1, 1, 4, 3, 2, 1)]
    dat = write_case(case_dir, "q4_invalid_detj", "Q4 invalid reversed node order", nodes, loads, bad_elements)
    expect_failure(exe, dat, "detJ <= 0")

    dat = case_dir / "q4_invalid_material.dat"
    lines = [
        "Q4 invalid material",
        "4 1 1 1",
        "1 1 1 1 0.0 0.0 0.0",
        "2 0 0 1 1.0 0.0 0.0",
        "3 0 0 1 1.0 1.0 0.0",
        "4 1 0 1 0.0 1.0 0.0",
        "1",
        "2",
        "2 1 5.0e-1",
        "3 1 5.0e-1",
        "2 1 1",
        "1 -1000.0 0.25 1.0 0",
        "1 1 2 3 4 1",
    ]
    dat.write_text("\n".join(lines) + "\n")
    expect_failure(exe, dat, "Young's modulus must be positive")

    (case_dir / "README.md").write_text(
        "# Q4 invalid-input checks\n\n"
        "These cases intentionally fail to verify clear diagnostics for active z DOF, reversed node order (`detJ <= 0`), and invalid material values.\n"
    )

def write_notes(results: Dict[str, object]) -> None:
    (ROOT / "data/q4_patch_single/README.md").write_text(
        "# Q4 single-element patch test\n\n"
        "Uniaxial constant-stress square patch. Expected stress: sigma_x=1, sigma_y=0, tau_xy=0.\n"
        f"Max stress error: {results['patch_single']['max_stress_error']:.6e}.\n"
        f"Max displacement error: {results['patch_single']['max_displacement_error']:.6e}.\n"
    )
    (ROOT / "data/q4_patch_multi/README.md").write_text(
        "# Q4 multi-element patch test\n\n"
        "2x2 constant-stress patch. Expected stress: sigma_x=1, sigma_y=0, tau_xy=0.\n"
        f"Max stress error: {results['patch_multi']['max_stress_error']:.6e}.\n"
        f"Max displacement error: {results['patch_multi']['max_displacement_error']:.6e}.\n"
    )
    (ROOT / "data/q4_patch_skew/README.md").write_text(
        "# Q4 skew-element patch test\n\n"
        "Single trapezoidal Q4 patch with nodes (0,0), (2,0), (2.5,1), (0,1). "
        "Expected stress: sigma_x=1, sigma_y=0, tau_xy=0.\n"
        f"Max stress error: {results['patch_skew']['max_stress_error']:.6e}.\n"
        f"Max displacement error: {results['patch_skew']['max_displacement_error']:.6e}.\n"
    )
    conv_lines = ["# Q4 cantilever convergence", "", "Reference is the 32x8 mesh tip displacement.", "", "| nx | ny | tip uy | rel error | observed order |", "|---:|---:|---:|---:|---:|"]
    for nx, ny, uy, err, order in results["convergence"]:
        order_s = "-" if math.isnan(order) else f"{order:.4f}"
        conv_lines.append(f"| {nx} | {ny} | {uy:.10e} | {err:.6e} | {order_s} |")
    conv_lines.extend(
        [
            "",
            "## ParaView post-processing",
            "",
            "The 32x8 reference mesh is used for the report displacement and von Mises contour figures. "
            "Generate the deformed-coordinate VTK file with a visual scale factor of 3:",
            "",
            "```bash",
            "python3 make/q4_to_vtk.py data/q4_convergence/q4_cantilever_32x8.dat data/q4_convergence/q4_cantilever_32x8.out data/q4_convergence/q4_cantilever_32x8.vtk --scale 3",
            "```",
        ]
    )
    (ROOT / "data/q4_convergence/README.md").write_text("\n".join(conv_lines) + "\n")
    (ROOT / "data/q4_validation/README.md").write_text(
        "# Q4 validation example\n\n"
        "The validation case reuses a 4x2 uniaxial membrane with an exact independent continuum solution. "
        "Expected stress is sigma_x=1, sigma_y=0, tau_xy=0 and displacement field is u=x/E, v=-nu*y/E.\n\n"
        f"Max stress error: {results['validation']['max_stress_error']:.6e}.\n"
        f"Max displacement error: {results['validation']['max_displacement_error']:.6e}.\n\n"
        "## ParaView post-processing\n\n"
        "`q4_validation_uniaxial.vtk` is a legacy VTK unstructured-grid file generated from the STAPpp `.dat/.out` pair and can be opened directly in ParaView. "
        "It contains deformed coordinates (scale 100), nodal displacement vectors, nodal `displacement_magnitude`, cell-average `sigma_x`, and cell-average plane-stress von Mises values.\n\n"
        "Regenerate with:\n\n"
        "```bash\n"
        "python3 make/q4_to_vtk.py data/q4_validation/q4_validation_uniaxial.dat data/q4_validation/q4_validation_uniaxial.out data/q4_validation/q4_validation_uniaxial.vtk --scale 100\n"
        "```\n"
    )


def validate_bar(exe: Path) -> None:
    import tempfile
    source = ROOT / "data/bar_examples/truss.dat"
    with tempfile.TemporaryDirectory(prefix="stappp-bar-") as tmp:
        dat = Path(tmp) / "truss.dat"
        dat.write_text(source.read_text())
        out = run_case(exe, dat)
        text = out.read_text()
    required = ["6.28230e-05", "1.66967e-03", "-4.88032e+08", "-3.45098e+08"]
    missing = [v for v in required if v not in text]
    if missing:
        raise AssertionError(f"bar regression missing expected values: {missing}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, type=Path)
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.exists():
        raise SystemExit(f"executable not found: {exe}")

    validate_bar(exe)
    results: Dict[str, object] = {}
    results["patch_single"] = check_patch(exe, ROOT / "data/q4_patch_single", "q4_patch_single", 1, 1, 1e-8, 1e-8)
    results["patch_multi"] = check_patch(exe, ROOT / "data/q4_patch_multi", "q4_patch_multi", 2, 2, 1e-6, 1e-8)
    results["patch_skew"] = check_skew_patch(exe)
    results["validation"] = check_patch(exe, ROOT / "data/q4_validation", "q4_validation_uniaxial", 4, 2, 1e-8, 1e-8)
    results["convergence"] = run_convergence(exe)
    run_invalid_cases(exe)
    write_notes(results)

    print("PASS: Bar regression")
    print(f"PASS: q4_patch_single max_stress_error={results['patch_single']['max_stress_error']:.6e} max_displacement_error={results['patch_single']['max_displacement_error']:.6e}")
    print(f"PASS: q4_patch_multi max_stress_error={results['patch_multi']['max_stress_error']:.6e} max_displacement_error={results['patch_multi']['max_displacement_error']:.6e}")
    print(f"PASS: q4_patch_skew max_stress_error={results['patch_skew']['max_stress_error']:.6e} max_displacement_error={results['patch_skew']['max_displacement_error']:.6e}")
    print(f"PASS: q4_validation max_stress_error={results['validation']['max_stress_error']:.6e} max_displacement_error={results['validation']['max_displacement_error']:.6e}")
    print("PASS: q4_convergence")
    for nx, ny, uy, err, order in results["convergence"]:
        order_s = "nan" if math.isnan(order) else f"{order:.6e}"
        print(f"  nx={nx} ny={ny} tip_uy={uy:.10e} rel_error={err:.6e} observed_order={order_s}")
    print("PASS: q4_invalid active_z/detJ/material diagnostics")
    return 0


if __name__ == "__main__":
    sys.exit(main())
