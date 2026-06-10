#!/usr/bin/env pvpython
"""Render Q4 VTK results with ParaView for the course-project report."""

from __future__ import annotations

import argparse
from pathlib import Path

from paraview.simple import (  # type: ignore
    ColorBy,
    CreateView,
    Delete,
    GetColorTransferFunction,
    GetScalarBar,
    HideScalarBarIfNotNeeded,
    LegacyVTKReader,
    LoadPalette,
    Render,
    ResetSession,
    SaveScreenshot,
    SetActiveView,
    Show,
)


ROOT = Path(__file__).resolve().parents[1]

JOBS = [
    {
        "vtk": ROOT / "data/q4_validation/q4_validation_uniaxial.vtk",
        "png": ROOT / "doc/figures/q4_validation_sigma_x.png",
        "association": "CELLS",
        "field": "sigma_x",
        "label": "sigma_x",
    },
    {
        "vtk": ROOT / "data/q4_patch_shear/q4_patch_shear.vtk",
        "png": ROOT / "doc/figures/q4_patch_shear_tau_xy.png",
        "association": "CELLS",
        "field": "tau_xy",
        "label": "tau_xy",
    },
    {
        "vtk": ROOT / "data/q4_plate_hole/q4_plate_hole_tension.vtk",
        "png": ROOT / "doc/figures/q4_plate_hole_von_mises.png",
        "association": "CELLS",
        "field": "von_mises_plane_stress",
        "label": "von Mises",
    },
    {
        "vtk": ROOT / "data/q4_convergence/q4_cantilever_32x8.vtk",
        "png": ROOT / "doc/figures/q4_cantilever_displacement.png",
        "association": "POINTS",
        "field": "displacement_magnitude",
        "label": "|u|",
    },
    {
        "vtk": ROOT / "data/q4_convergence/q4_cantilever_32x8.vtk",
        "png": ROOT / "doc/figures/q4_cantilever_von_mises.png",
        "association": "CELLS",
        "field": "von_mises_plane_stress",
        "label": "von Mises",
    },
]


def set_xy_camera(view, bounds: tuple[float, float, float, float, float, float], width: int, height: int) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
    span_x = max(xmax - xmin, 1.0e-12)
    span_y = max(ymax - ymin, 1.0e-12)
    aspect = width / height

    view.CameraParallelProjection = 1
    view.CameraPosition = [cx, cy, cz + max(span_x, span_y, 1.0) * 4.0]
    view.CameraFocalPoint = [cx, cy, cz]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = 0.58 * max(span_y, span_x / aspect)


def render_job(job: dict[str, object], width: int, height: int) -> None:
    vtk = Path(job["vtk"])
    png = Path(job["png"])
    field = str(job["field"])
    association = str(job["association"])
    label = str(job["label"])

    if not vtk.exists():
        raise FileNotFoundError(vtk)
    png.parent.mkdir(parents=True, exist_ok=True)

    ResetSession()
    LoadPalette("WhiteBackground")
    view = CreateView("RenderView")
    SetActiveView(view)
    view.ViewSize = [width, height]
    try:
        view.UseColorPaletteForBackground = 0
        view.BackgroundColorMode = "Single Color"
    except Exception:
        pass
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 1

    reader = LegacyVTKReader(registrationName=vtk.stem, FileNames=[str(vtk)])
    reader.UpdatePipeline()

    display = Show(reader, view, "UnstructuredGridRepresentation")
    display.Representation = "Surface With Edges"
    display.EdgeColor = [0.08, 0.08, 0.08]
    display.LineWidth = 1.0

    ColorBy(display, (association, field))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)

    lut = GetColorTransferFunction(field)
    try:
        lut.ApplyPreset("Viridis (matplotlib)", True)
    except Exception:
        lut.ApplyPreset("Cool to Warm", True)

    scalar_bar = GetScalarBar(lut, view)
    scalar_bar.Title = label
    scalar_bar.ComponentTitle = ""
    scalar_bar.TitleFontSize = 10
    scalar_bar.LabelFontSize = 9
    scalar_bar.DrawTickMarks = 1

    set_xy_camera(view, reader.GetDataInformation().GetBounds(), width, height)
    Render(view)
    SaveScreenshot(str(png), view, ImageResolution=[width, height], TransparentBackground=0)
    HideScalarBarIfNotNeeded(lut, view)
    Delete(reader)
    Delete(view)
    print(f"wrote {png.relative_to(ROOT)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=1100)
    args = parser.parse_args()

    for job in JOBS:
        render_job(job, args.width, args.height)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
