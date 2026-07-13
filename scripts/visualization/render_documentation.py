#!/usr/bin/env python3
"""Render AdditiveFOAM documentation figures with ParaView's Python interface.

Run this script with ``pvbatch`` after the tutorial cases have been reconstructed.
Each case directory must contain an empty ``case.foam`` marker file.
"""

from pathlib import Path
import argparse
import shutil
import subprocess

from PIL import Image
from paraview.simple import (  # type: ignore
    _DisableFirstRenderCameraReset,
    CellDatatoPointData,
    ColorBy,
    Contour,
    CreateLayout,
    CreateView,
    Delete,
    GetColorTransferFunction,
    GetScalarBar,
    OpenFOAMReader,
    Render,
    ResetSession,
    SaveScreenshot,
    SetActiveView,
    Show,
    Slice,
    Text,
    Threshold,
)


_DisableFirstRenderCameraReset()

IMAGE_WIDTH = 1600
ORNL_GREEN = [0.0, 0.475, 0.325]
DARK_GREEN = [0.0, 0.235, 0.165]
TEXT_COLOR = [0.12, 0.12, 0.12]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-case", type=Path, required=True)
    parser.add_argument("--multi-beam-case", type=Path, required=True)
    parser.add_argument("--nlight-case", type=Path, required=True)
    parser.add_argument("--tabulated-case", type=Path, required=True)
    parser.add_argument("--amr-case", type=Path, required=True)
    parser.add_argument("--multi-layer-case", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--only",
        choices=(
            "all",
            "amb2018",
            "multi-beam",
            "nlight-afx",
            "tabulated",
            "multi-layer",
            "amr",
            "multi-layer-fields",
        ),
        default="all",
        help="Render one asset group in a fresh ParaView process (default: all).",
    )
    return parser.parse_args()


def reader(case, arrays):
    source = OpenFOAMReader(FileName=str(case / "case.foam"))
    source.MeshRegions = ["internalMesh"]
    source.CellArrays = arrays
    source.UpdatePipelineInformation()
    return source


def render_view(height=900):
    view = CreateView("RenderView")
    view.ViewSize = [IMAGE_WIDTH, height]
    view.Background = [1.0, 1.0, 1.0]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    view.AxesGrid.Visibility = 0
    if hasattr(view, "UseDistributedRenderingForRender"):
        view.UseDistributedRenderingForRender = 0
    if hasattr(view, "UseDistributedRenderingForLODRender"):
        view.UseDistributedRenderingForLODRender = 0
    return view


def set_isometric_camera(view):
    view.CameraPosition = [0.0036, -0.0032, 0.0023]
    view.CameraFocalPoint = [0.0010, 0.0, -0.00010]
    view.CameraViewUp = [-0.35, 0.31, 0.88]
    view.CameraParallelScale = 0.00145


def set_top_camera(view, parallel_scale=0.00060):
    view.CameraPosition = [0.0010, 0.0, 0.01]
    view.CameraFocalPoint = [0.0010, 0.0, -0.00005]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = parallel_scale


def add_text(view, value, position="Upper Left Corner", size=32):
    text = Text(Text=value)
    display = Show(text, view)
    display.WindowLocation = position
    display.FontSize = size
    display.Color = TEXT_COLOR
    return text, display


def style_scalar_bar(bar, length=0.42, font_size=26):
    bar.Orientation = "Horizontal"
    bar.WindowLocation = "Lower Center"
    bar.TitleColor = TEXT_COLOR
    bar.LabelColor = TEXT_COLOR
    bar.TitleFontSize = font_size
    bar.LabelFontSize = font_size - 2
    bar.ScalarBarLength = length
    bar.ScalarBarThickness = 24


def temperature_lut(
    view, display, maximum, association="CELLS", custom_labels=None
):
    ColorBy(display, (association, "T"))
    lut = GetColorTransferFunction("T")
    lut.RGBPoints = [
        300.0, 0.230, 0.299, 0.754,
        0.25 * (maximum - 300.0) + 300.0, 0.554, 0.690, 0.996,
        0.50 * (maximum - 300.0) + 300.0, 0.865, 0.865, 0.865,
        0.75 * (maximum - 300.0) + 300.0, 0.957, 0.598, 0.478,
        maximum, 0.706, 0.016, 0.150,
    ]
    lut.ColorSpace = "RGB"
    lut.RescaleTransferFunction(300.0, maximum)
    bar = GetScalarBar(lut, view)
    bar.Title = "Temperature"
    bar.ComponentTitle = "K"
    bar.UseCustomLabels = 1
    bar.CustomLabels = custom_labels or [300.0, 1000.0, 1620.0, maximum]
    style_scalar_bar(bar, length=0.44, font_size=26)
    display.SetScalarBarVisibility(view, True)
    return lut


def encode_animation(frames, frame_prefix, poster_index, poster, video):
    shutil.copyfile(
        frames / "{}-{:04d}.png".format(frame_prefix, poster_index), poster
    )
    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-loglevel",
            "error",
            "-framerate",
            "1",
            "-i",
            str(frames / "{}-%04d.png".format(frame_prefix)),
            "-vf",
            "fps=30,format=yuv420p",
            "-c:v",
            "libx264",
            "-preset",
            "slow",
            "-tune",
            "stillimage",
            "-crf",
            "10",
            "-profile:v",
            "high",
            "-level",
            "4.1",
            "-movflags",
            "+faststart",
            str(video),
        ],
        check=True,
    )
    shutil.rmtree(frames)


def render_temperature_animation(
    case,
    output_dir,
    slug,
    title,
    alloy,
    solidus_temperature,
    liquidus_temperature,
    maximum_temperature,
    poster_time,
    poster_name=None,
):
    source = reader(case, ["T", "qDot"])
    times = list(source.TimestepValues)
    source.UpdatePipeline(times[0])
    bounds = source.GetDataInformation().GetBounds()
    # Place the cutter infinitesimally inside the mesh so vtkCutter samples
    # the exposed boundary without losing coplanar top faces numerically.
    slice_z = bounds[5] - 1.0e-6 * (bounds[5] - bounds[4])
    frames = output_dir / "{}-frames".format(slug)
    frames.mkdir(parents=True, exist_ok=True)

    point_data = CellDatatoPointData(Input=source)
    point_data.PassCellData = 1
    section = Slice(Input=point_data)
    section.SliceType = "Plane"
    section.SliceType.Origin = [0.0010, 0.0, slice_z]
    section.SliceType.Normal = [0.0, 0.0, 1.0]
    solidus = Contour(Input=section)
    solidus.ContourBy = ["POINTS", "T"]
    solidus.Isosurfaces = [solidus_temperature]
    liquidus = Contour(Input=section)
    liquidus.ContourBy = ["POINTS", "T"]
    liquidus.Isosurfaces = [liquidus_temperature]

    view = render_view(height=600)
    display = Show(section, view)
    display.Representation = "Surface"
    display.Interpolation = "Gouraud"
    temperature_lut(
        view,
        display,
        maximum=maximum_temperature,
        association="POINTS",
        custom_labels=[300.0, maximum_temperature],
    )
    solidus_display = Show(solidus, view)
    solidus_display.Representation = "Surface"
    solidus_display.DiffuseColor = [0.0, 0.0, 0.0]
    solidus_display.AmbientColor = [0.0, 0.0, 0.0]
    solidus_display.LineWidth = 3.0
    solidus_display.ColorArrayName = [None, ""]
    liquidus_display = Show(liquidus, view)
    liquidus_display.Representation = "Surface"
    liquidus_display.DiffuseColor = [1.0, 1.0, 1.0]
    liquidus_display.AmbientColor = [1.0, 1.0, 1.0]
    liquidus_display.LineWidth = 3.0
    liquidus_display.ColorArrayName = [None, ""]
    header, header_display = add_text(view, "", size=28)
    contour_key, contour_key_display = add_text(
        view,
        "black: {:.0f} K solidus\nwhite: {:.0f} K liquidus".format(
            solidus_temperature, liquidus_temperature
        ),
        position="Upper Right Corner",
        size=24,
    )
    set_top_camera(view, parallel_scale=0.00060)
    SetActiveView(view)

    for index, time in enumerate(times):
        source.UpdatePipeline(time)
        bounds = source.GetDataInformation().GetBounds()
        slice_z = bounds[5] - 1.0e-6 * (bounds[5] - bounds[4])
        section.SliceType.Origin = [0.0010, 0.0, slice_z]
        section.UpdatePipeline(time)
        solidus.UpdatePipeline(time)
        liquidus.UpdatePipeline(time)
        view.ViewTime = time
        header.Text = "{}  |  {}\nTop-surface T(x,y)  |  t = {:.1f} ms".format(
            title, alloy, time * 1.0e3
        )
        # ParaView can reapply its global palette as views are activated. Reset
        # the canvas explicitly so every frame has the same white background.
        view.UseColorPaletteForBackground = 0
        view.Background = [1.0, 1.0, 1.0]
        Render(view)
        SaveScreenshot(
            str(frames / "{}-{:04d}.png".format(slug, index)),
            view,
            OverrideColorPalette="WhiteBackground",
            TransparentBackground=0,
        )

    poster_index = min(
        range(len(times)), key=lambda index: abs(times[index] - poster_time)
    )
    encode_animation(
        frames,
        slug,
        poster_index,
        output_dir / (poster_name or "{}.png".format(slug)),
        output_dir / "{}.mp4".format(slug),
    )
    # Each animation owns its view and pipeline. Releasing them prevents
    # off-screen framebuffer tiles from leaking into the next render.
    for proxy in (
        contour_key_display,
        contour_key,
        header_display,
        header,
        liquidus_display,
        liquidus,
        solidus_display,
        solidus,
        display,
        section,
        point_data,
        source,
        view,
    ):
        Delete(proxy)
    ResetSession()


def render_amr(case, output_dir):
    frames = output_dir / "amr-frames"
    frames.mkdir(parents=True, exist_ok=True)
    source = reader(case, ["cellLevel", "refinementField"])

    view = render_view(height=600)
    mesh = Show(source, view)
    mesh.Representation = "Surface With Edges"
    mesh.EdgeColor = [0.30, 0.30, 0.30]
    mesh.LineWidth = 0.6
    ColorBy(mesh, ("CELLS", "cellLevel"))
    level_lut = GetColorTransferFunction("cellLevel")
    level_lut.RGBPoints = [
        0.0, 0.92, 0.92, 0.92,
        1.0, DARK_GREEN[0], DARK_GREEN[1], DARK_GREEN[2],
    ]
    level_lut.ColorSpace = "RGB"
    level_lut.RescaleTransferFunction(0.0, 1.0)
    level_bar = GetScalarBar(level_lut, view)
    level_bar.Title = "Cell refinement level"
    level_bar.ComponentTitle = ""
    level_bar.UseCustomLabels = 1
    level_bar.CustomLabels = [0.0, 1.0]
    style_scalar_bar(level_bar, length=0.30, font_size=24)
    mesh.SetScalarBarVisibility(view, True)

    requested = Threshold(Input=source)
    requested.Scalars = ["CELLS", "refinementField"]
    requested.LowerThreshold = 0.5
    requested.UpperThreshold = 1.5
    requested.ThresholdMethod = "Between"
    requested_display = Show(requested, view)
    requested_display.Representation = "Surface With Edges"
    requested_display.DiffuseColor = ORNL_GREEN
    requested_display.AmbientColor = ORNL_GREEN
    requested_display.EdgeColor = ORNL_GREEN
    requested_display.Opacity = 0.45
    requested_display.LineWidth = 1.2
    requested_display.ColorArrayName = [None, ""]

    title, title_display = add_text(
        view,
        "Dynamic AMR  |  IN625\nGreen: refinementField = 1",
        size=28,
    )
    time_text, time_display = add_text(
        view,
        "t = 0.0 ms",
        position="Upper Right Corner",
        size=28,
    )
    set_top_camera(view)
    SetActiveView(view)

    times = list(source.TimestepValues)
    for index, time in enumerate(times):
        source.UpdatePipeline(time)
        requested.UpdatePipeline(time)
        view.ViewTime = time
        time_text.Text = "t = {:.1f} ms".format(time * 1000.0)
        Render(view)
        SaveScreenshot(
            str(frames / "amr-{:04d}.png".format(index)),
            view,
            OverrideColorPalette="WhiteBackground",
            TransparentBackground=0,
        )

    encode_animation(
        frames,
        "amr",
        1,
        output_dir / "amr-refinement.png",
        output_dir / "amr-refinement.mp4",
    )


def configure_panel(view, source, field, title):
    display = Show(source, view)
    display.Representation = "Surface"
    display.Interpolation = "Gouraud"
    if field == "T":
        temperature_lut(view, display, maximum=3300.0)
    else:
        ColorBy(display, ("CELLS", field))
        lut = GetColorTransferFunction(field)
        lut.RGBPoints = [
            0.0, 0.02, 0.02, 0.02,
            1.0, 0.94, 0.94, 0.94,
        ]
        lut.ColorSpace = "RGB"
        lut.RescaleTransferFunction(0.0, 1.0)
        bar = GetScalarBar(lut, view)
        bar.Title = "Powder fraction"
        bar.ComponentTitle = ""
        bar.UseCustomLabels = 1
        bar.CustomLabels = [0.0, 0.5, 1.0]
        style_scalar_bar(bar, length=0.38, font_size=24)
        display.SetScalarBarVisibility(view, True)
    add_text(view, title, size=32)
    set_top_camera(view, parallel_scale=0.00045)
    return display


def render_multi_layer(case, output_dir):
    source = reader(case, ["T", "alpha.powder", "alpha.solid"])
    time = 0.0225
    source.UpdatePipeline(time)

    top = render_view(height=450)
    bottom = render_view(height=450)
    top.ViewTime = time
    bottom.ViewTime = time
    configure_panel(top, source, "T", "Temperature  |  IN625  |  t = 22.5 ms")
    configure_panel(
        bottom,
        source,
        "alpha.powder",
        "Powder fraction  |  IN625  |  layer 2",
    )

    layout = CreateLayout(name="Multi-layer fields")
    layout.SplitVertical(0, 0.5)
    layout.AssignView(1, top)
    layout.AssignView(2, bottom)
    layout.SetSize(IMAGE_WIDTH, 900)
    Render(top)
    Render(bottom)
    SaveScreenshot(
        str(output_dir / "multi-layer-fields.png"),
        layout,
        OverrideColorPalette="WhiteBackground",
        TransparentBackground=0,
    )


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.only in ("all", "amb2018"):
        render_temperature_animation(
            args.baseline_case.resolve(),
            args.output_dir.resolve(),
            "amb2018-temperature",
            "AMB2018-02-B",
            "IN625",
            1410.0,
            1620.0,
            2700.0,
            0.002,
            poster_name="quick-start-temperature.png",
        )
    if args.only in ("all", "multi-beam"):
        render_temperature_animation(
            args.multi_beam_case.resolve(),
            args.output_dir.resolve(),
            "multi-beam-temperature",
            "Multi-Beam",
            "IN625",
            1410.0,
            1620.0,
            3100.0,
            0.002,
        )
    if args.only in ("all", "nlight-afx"):
        render_temperature_animation(
            args.nlight_case.resolve(),
            args.output_dir.resolve(),
            "nlight-afx-temperature",
            "nLight AFX",
            "SS316L",
            1471.0,
            1709.0,
            2900.0,
            0.002,
        )
    if args.only in ("all", "tabulated"):
        render_temperature_animation(
            args.tabulated_case.resolve(),
            args.output_dir.resolve(),
            "tabulated-temperature",
            "Tabulated Beam Profile",
            "AlSi10Mg",
            850.0,
            870.0,
            1600.0,
            0.002,
        )
    if args.only in ("all", "multi-layer"):
        render_temperature_animation(
            args.multi_layer_case.resolve(),
            args.output_dir.resolve(),
            "multi-layer-temperature",
            "Multi-Layer PBF",
            "IN625",
            1410.0,
            1620.0,
            3300.0,
            0.0225,
        )
    if args.only in ("all", "amr"):
        render_amr(args.amr_case.resolve(), args.output_dir.resolve())
    if args.only in ("all", "multi-layer-fields"):
        render_multi_layer(
            args.multi_layer_case.resolve(), args.output_dir.resolve()
        )
    for figure in args.output_dir.glob("*.png"):
        with Image.open(figure) as image:
            image.save(figure, dpi=(300, 300))


if __name__ == "__main__":
    main()
