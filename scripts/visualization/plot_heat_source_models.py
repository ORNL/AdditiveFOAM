#!/usr/bin/env python3
"""Plot normalized AdditiveFOAM heat-source distributions."""

from pathlib import Path
import argparse
from math import erf, pi, sqrt

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from plot_style import PUBLICATION_STYLE


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tabulated-profile", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--projection-output", type=Path, required=True)
    return parser.parse_args()


def read_tabulated(path):
    with path.open(encoding="utf-8") as stream:
        nx, ny = (int(value) for value in stream.readline().split())
        x0, y0 = (float(value) for value in stream.readline().split())
        dx, dy = (float(value) for value in stream.readline().split())
        values = np.loadtxt(stream)
    values = np.asarray(values).reshape(ny, nx)
    x = x0 + dx * np.arange(nx)
    y = y0 + dy * np.arange(ny)
    return x, y, values


def ring_intensity(radius, sigma, radial_coordinate):
    return np.exp(-0.5 * ((radial_coordinate - radius) / sigma) ** 2) + np.exp(
        -0.5 * ((radial_coordinate + radius) / sigma) ** 2
    )


def ring_area(radius, sigma):
    return 2.0 * pi * sigma * (
        2.0 * sigma * np.exp(-(radius**2) / (2.0 * sigma**2))
        + sqrt(2.0 * pi) * radius * erf(radius / (sqrt(2.0) * sigma))
    )


def normalize(top, section):
    maximum = max(float(np.nanmax(top)), float(np.nanmax(section)))
    return top / maximum, section / maximum


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    lateral = np.linspace(-1.5, 1.5, 401)
    depth = np.linspace(0.0, 1.5, 241)
    xx, yy = np.meshgrid(lateral, lateral)
    xs, zz = np.meshgrid(lateral, depth)

    models = []

    exponent = 4.0
    scale = 1.0 / 2.0 ** (1.0 / exponent)
    top = np.exp(-((xx / scale) ** 2 + (yy / scale) ** 2) ** (exponent / 2.0))
    section = np.exp(
        -((xs / scale) ** 2 + (zz / scale) ** 2) ** (exponent / 2.0)
    )
    models.append(("superGaussian\n$k=4$", *normalize(top, section)))

    exponent = 7.95
    contraction = 2.72
    scale = 1.0 / 2.0 ** (1.0 / exponent)
    top = np.exp(-((xx / scale) ** 2 + (yy / scale) ** 2) ** (exponent / 2.0))
    g = np.zeros_like(zz)
    inside = zz < 1.0
    g[inside] = (1.0 - zz[inside] ** contraction) ** (1.0 / contraction)
    section = np.zeros_like(zz)
    section[inside] = np.exp(
        -((xs[inside] / (scale * g[inside])) ** 2) ** (exponent / 2.0)
    )
    models.append(
        ("modifiedSuperGaussian\n$k=7.95$, $m=2.72$", *normalize(top, section))
    )

    projection = np.exp(-3.0 * zz**2)
    top = np.exp(-2.0 * (xx**2 + yy**2))
    section = np.exp(-2.0 * xs**2) * projection
    models.append(("projectedGaussian\n$k_z=2$", *normalize(top, section)))

    # Tutorial Index 6 parameters, normalized by its lateral dimension.
    alpha = 0.902
    dimension = 1.42905e-4
    inner_radius, inner_sigma = 3.3500e-5 / dimension, 1.8510e-5 / dimension
    outer_radius, outer_sigma = 1.0141e-4 / dimension, 1.63475e-5 / dimension
    radial_top = np.sqrt(xx**2 + yy**2)
    radial_section = np.abs(xs)
    top = (1.0 - alpha) * ring_intensity(
        inner_radius, inner_sigma, radial_top
    ) / ring_area(inner_radius, inner_sigma)
    top += alpha * ring_intensity(
        outer_radius, outer_sigma, radial_top
    ) / ring_area(outer_radius, outer_sigma)
    lateral_section = (1.0 - alpha) * ring_intensity(
        inner_radius, inner_sigma, radial_section
    ) / ring_area(inner_radius, inner_sigma)
    lateral_section += alpha * ring_intensity(
        outer_radius, outer_sigma, radial_section
    ) / ring_area(outer_radius, outer_sigma)
    section = lateral_section * projection
    models.append(("nLightAFX Index 6\n$\\alpha=0.902$", *normalize(top, section)))

    table_x, table_y, table = read_tabulated(args.tabulated_profile)
    table_dimension = 2.50e-4
    table_x = table_x / table_dimension
    table_y = table_y / table_dimension
    row = table[np.argmin(np.abs(table_y))]
    lateral_table = np.interp(lateral, table_x, row, left=0.0, right=0.0)
    section = lateral_table[np.newaxis, :] * np.exp(-3.0 * depth[:, np.newaxis] ** 2)
    models.append(("tabulated Index 6\nmeasured profile", *normalize(table, section)))

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    figure = plt.figure(figsize=(8.5, 15.2), facecolor="white")
    grid = figure.add_gridspec(
        5,
        2,
        left=0.32,
        right=0.97,
        bottom=0.07,
        top=0.95,
        width_ratios=(1.0, 1.0),
        hspace=0.04,
        wspace=0.12,
    )
    row_labels = (
        "(a)  superGaussian\n      $k=4$",
        "(b)  modifiedSuperGaussian\n      $k=7.95$, $m=2.72$",
        "(c)  projectedGaussian\n      $k_z=2$",
        "(d)  nLightAFX Index 6\n      $\\alpha=0.902$",
        "(e)  tabulated Index 6\n      measured profile",
    )
    images = []
    panel_axes = []

    for row_index, ((label, top, section), row_label) in enumerate(
        zip(models, row_labels)
    ):
        top_axis = figure.add_subplot(grid[row_index, 0])
        section_axis = figure.add_subplot(grid[row_index, 1])
        top_extent = [-1.5, 1.5, -1.5, 1.5]
        if label.startswith("tabulated"):
            top_extent = [table_x[0], table_x[-1], table_y[0], table_y[-1]]
        images.append(
            top_axis.imshow(
                top,
                origin="lower",
                extent=top_extent,
                vmin=0.0,
                vmax=1.0,
                cmap="viridis",
                interpolation="bilinear",
                aspect="equal",
            )
        )
        section_axis.imshow(
            section,
            origin="upper",
            extent=[-1.5, 1.5, 1.5, 0.0],
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
            interpolation="bilinear",
            aspect="auto",
        )
        top_axis.set_xlim(-1.5, 1.5)
        top_axis.set_ylim(-1.5, 1.5)
        section_axis.set_xlim(-1.5, 1.5)
        section_axis.set_ylim(1.5, 0.0)
        top_axis.set_box_aspect(1.0)
        section_axis.set_box_aspect(1.0)
        top_axis.set_xticks([])
        top_axis.set_yticks([])
        section_axis.set_xticks([])
        section_axis.set_yticks([])
        panel_axes.append((top_axis, section_axis))
        centre = top_axis.get_position().y0 + top_axis.get_position().height / 2
        figure.text(
            0.025,
            centre,
            row_label,
            ha="left",
            va="center",
            fontsize=10.5,
            linespacing=1.35,
        )

    projected_top = max(
        axis.get_position().y1 for axis in panel_axes[2]
    )
    projected_bottom = min(
        axis.get_position().y0 for axis in panel_axes[4]
    )
    projected_outline = Rectangle(
        (0.012, projected_bottom - 0.006),
        0.976,
        projected_top - projected_bottom + 0.010,
        transform=figure.transFigure,
        fill=False,
        edgecolor="#007A53",
        linewidth=1.4,
        linestyle="--",
        clip_on=False,
    )
    figure.add_artist(projected_outline)
    figure.text(
        0.012,
        0.5 * (projected_top + projected_bottom),
        "Common projected-source formulation",
        ha="center",
        va="center",
        rotation=90,
        fontsize=10.5,
        color="#006B49",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": "white",
            "edgecolor": "#007A53",
            "linewidth": 1.0,
        },
    )

    first_top = figure.axes[0].get_position()
    first_section = figure.axes[1].get_position()
    figure.text(
        first_top.x0 + first_top.width / 2,
        0.965,
        "Top surface",
        ha="center",
        va="center",
        fontsize=12,
    )
    figure.text(
        first_section.x0 + first_section.width / 2,
        0.965,
        "Centre-plane section",
        ha="center",
        va="center",
        fontsize=12,
    )
    colorbar_axis = figure.add_axes([0.40, 0.025, 0.45, 0.012])
    colorbar = figure.colorbar(images[0], cax=colorbar_axis, orientation="horizontal")
    colorbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    colorbar.set_label("Normalized source magnitude", labelpad=2)
    figure.savefig(args.output, dpi=300, facecolor="white")

    with plt.rc_context(PUBLICATION_STYLE):
        projection_figure, projection_axis = plt.subplots(
            figsize=(4.5, 3.0), constrained_layout=True
        )
        colors = ("#0072B2", "#E69F00", "#009E73", "#CC3311")
        styles = ("-", "--", "-.", ":")
        for exponent, color, style in zip((1, 2, 4, 8), colors, styles):
            projection_axis.plot(
                depth,
                np.exp(-3.0 * depth**exponent),
                color=color,
                linestyle=style,
                label=f"$k={exponent}$",
            )
        projection_axis.set(
            xlabel="Normalized depth,  $z/d_z$",
            ylabel="$p(z)$",
            xlim=(0.0, 1.5),
            ylim=(0.0, 1.03),
        )
        projection_axis.grid(True, which="both", linestyle="-", alpha=0.1)
        projection_axis.minorticks_on()
        projection_axis.legend(ncol=4, frameon=True, loc="upper right")
        projection_figure.savefig(args.projection_output)
        plt.close(projection_figure)


if __name__ == "__main__":
    main()
