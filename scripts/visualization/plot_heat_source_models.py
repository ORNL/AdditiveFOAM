#!/usr/bin/env python3
"""Plot normalized AdditiveFOAM heat-source distributions."""

from pathlib import Path
import argparse
from math import erf, pi, sqrt
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Rectangle

from plot_style import PUBLICATION_STYLE


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tabulated-profile", type=Path, required=True)
    parser.add_argument("--nlight-config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--projection-output", type=Path, required=True)
    parser.add_argument("--metrics-output", type=Path, required=True)
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


def read_nlight_mode(path, mode="Index6"):
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"\b{re.escape(mode)}\s*\{{(.*?)\}}", text, re.DOTALL)
    if match is None:
        raise ValueError(f"Could not find {mode} in {path}")
    block = match.group(1)
    values = {}
    for name in ("alpha", "r0", "sigma0", "r1", "sigma1"):
        entry = re.search(rf"\b{name}\s+([^;]+);", block)
        if entry is None:
            raise ValueError(f"Could not find {name} in {mode} of {path}")
        values[name] = float(entry.group(1))
    return values


def tabulated_metrics(x, y, values):
    """Exact raw moments of the nodal table's bilinear interpolant."""
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    moments = np.zeros(6)
    for j in range(len(y) - 1):
        for i in range(len(x) - 1):
            f00, f10 = values[j, i], values[j, i + 1]
            f01, f11 = values[j + 1, i], values[j + 1, i + 1]
            m00 = 0.25 * (f00 + f10 + f01 + f11)
            m10 = (f00 + 2 * f10 + f01 + 2 * f11) / 12
            m01 = (f00 + f10 + 2 * f01 + 2 * f11) / 12
            m20 = (f00 + 3 * f10 + f01 + 3 * f11) / 24
            m02 = (f00 + f10 + 3 * f01 + 3 * f11) / 24
            m11 = (f00 + 2 * f10 + 2 * f01 + 4 * f11) / 36
            x0, y0 = x[i], y[j]
            da = dx * dy
            moments += da * np.array(
                [
                    m00,
                    x0 * m00 + dx * m10,
                    y0 * m00 + dy * m01,
                    x0**2 * m00 + 2 * x0 * dx * m10 + dx**2 * m20,
                    y0**2 * m00 + 2 * y0 * dy * m01 + dy**2 * m02,
                    x0 * y0 * m00
                    + x0 * dy * m01
                    + y0 * dx * m10
                    + dx * dy * m11,
                ]
            )
    m00, m10, m01, m20, m02, m11 = moments
    centroid = np.array([m10 / m00, m01 / m00])
    variance_x = m20 / m00 - centroid[0] ** 2
    variance_y = m02 / m00 - centroid[1] ** 2
    covariance = m11 / m00 - centroid[0] * centroid[1]
    eigenvalues = np.linalg.eigvalsh(
        np.array([[variance_x, covariance], [covariance, variance_y]])
    )
    major, minor = 4 * np.sqrt(eigenvalues[::-1])
    return centroid, sqrt(major * minor)


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


def plot_profile_metrics(output):
    """Connect principal beam metrics to the shared reference aspect ratio."""
    with plt.rc_context(PUBLICATION_STYLE):
        figure, (profile_axis, flow_axis) = plt.subplots(
            1, 2, figsize=(8.0, 3.8), constrained_layout=True
        )

        centroid = np.array([0.25, -0.18])
        sigma_major = 0.72
        sigma_minor = 0.36
        theta = np.deg2rad(28.0)
        major_axis = np.array([np.cos(theta), np.sin(theta)])
        minor_axis = np.array([-np.sin(theta), np.cos(theta)])
        rotation = np.column_stack((major_axis, minor_axis))
        covariance = rotation @ np.diag([sigma_major**2, sigma_minor**2]) @ rotation.T

        x = np.linspace(-1.5, 2.0, 360)
        y = np.linspace(-1.5, 1.5, 320)
        xx, yy = np.meshgrid(x, y)
        offsets = np.stack((xx - centroid[0], yy - centroid[1]), axis=-1)
        quadratic = np.einsum(
            "...i,ij,...j->...", offsets, np.linalg.inv(covariance), offsets
        )
        profile_axis.contourf(
            xx, yy, np.exp(-0.5 * quadratic), levels=24, cmap="viridis"
        )
        profile_axis.add_patch(
            Ellipse(
                centroid,
                width=4.0 * sigma_major,
                height=4.0 * sigma_minor,
                angle=np.rad2deg(theta),
                fill=False,
                edgecolor="white",
                linewidth=1.0,
                linestyle="--",
            )
        )
        for sigma, direction, color, label, offset in (
            (sigma_major, major_axis, "white", "$D_{4\\sigma,\\,major}$", 0.10),
            (sigma_minor, minor_axis, "#F0E442", "$D_{4\\sigma,\\,minor}$", -0.14),
        ):
            endpoints = np.vstack(
                (centroid - 2.0 * sigma * direction, centroid + 2.0 * sigma * direction)
            )
            profile_axis.plot(endpoints[:, 0], endpoints[:, 1], color=color, linewidth=1.8)
            label_point = centroid + (1.15 * sigma + offset) * direction
            profile_axis.text(
                label_point[0],
                label_point[1],
                label,
                color=color,
                ha="center",
                va="bottom",
                fontsize=9,
                bbox={"facecolor": "black", "alpha": 0.25, "edgecolor": "none"},
            )
        profile_axis.scatter(*centroid, color="white", edgecolor="black", s=28, zorder=4)
        profile_axis.annotate(
            "centroid",
            xy=centroid,
            xytext=(centroid[0] - 0.72, centroid[1] - 0.55),
            color="white",
            arrowprops={"arrowstyle": "->", "color": "white"},
        )
        area_equivalent = np.sqrt((4.0 * sigma_major) * (4.0 * sigma_minor))
        profile_axis.text(
            0.02,
            0.98,
            "$D_{ref}=\\sqrt{D_{major}D_{minor}}$\n"
            f"area-equivalent width = {area_equivalent:.2f} (scaled)",
            transform=profile_axis.transAxes,
            ha="left",
            va="top",
            color="white",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "black", "alpha": 0.48},
        )
        profile_axis.set(
            title="Rotated elliptical profile metrics",
            xlabel="Beam-plane x",
            ylabel="Beam-plane y",
            xlim=(x.min(), x.max()),
            ylim=(y.min(), y.max()),
        )
        profile_axis.set_aspect("equal")

        flow_axis.axis("off")

        def box(x_pos, y_pos, text, color="#E8F1F5"):
            flow_axis.text(
                x_pos,
                y_pos,
                text,
                transform=flow_axis.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                bbox={
                    "boxstyle": "round,pad=0.38",
                    "facecolor": color,
                    "edgecolor": "0.25",
                },
            )

        def arrow(start, end):
            flow_axis.annotate(
                "",
                xy=end,
                xytext=start,
                xycoords=flow_axis.transAxes,
                textcoords=flow_axis.transAxes,
                arrowprops={"arrowstyle": "->", "color": "0.3", "linewidth": 1.1},
            )

        box(0.25, 0.88, "profileMetrics\n$D_{major}, D_{minor}$")
        box(0.25, 0.67, "D4Sigma selector\n$D_{ref}$", "#DDF2E9")
        box(0.75, 0.88, "depthReference\nconstant or isotherm")
        box(0.75, 0.67, "reference depth\n$d_{ref}$", "#DDF2E9")
        box(0.32, 0.43, "$a=2d_{ref}/D_{ref}$", "#FFF0CC")
        box(0.77, 0.43, "$d_{source}=\\max($\n$d_{configured},d_{ref})$\napplied source depth")
        box(0.25, 0.16, "exponential projection\n$n(a)$")
        box(0.75, 0.16, "Kelly absorption\n$\\eta(a)$")
        arrow((0.25, 0.82), (0.25, 0.73))
        arrow((0.75, 0.82), (0.75, 0.73))
        arrow((0.25, 0.62), (0.30, 0.48))
        arrow((0.70, 0.62), (0.39, 0.47))
        arrow((0.75, 0.62), (0.76, 0.50))
        arrow((0.30, 0.38), (0.26, 0.22))
        arrow((0.39, 0.39), (0.69, 0.20))
        flow_axis.set_title("Reference-dimension data flow")

        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output)
        plt.close(figure)


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
    models.append(
        (
            "projected\nsuperGaussian profile + exponential projection",
            *normalize(top, section),
        )
    )

    # Read the same shared Index 6 coefficients used by the tutorial.
    mode = read_nlight_mode(args.nlight_config)
    alpha = mode["alpha"]
    radial_grid = np.linspace(0.0, 4.0e-4, 20001)
    component0 = ring_intensity(mode["r0"], mode["sigma0"], radial_grid)
    component1 = ring_intensity(mode["r1"], mode["sigma1"], radial_grid)
    j0 = np.trapezoid(component0 * radial_grid, radial_grid)
    j1 = np.trapezoid(component1 * radial_grid, radial_grid)
    combined = (1.0 - alpha) * component0 + alpha * component1 * j0 / j1
    radial_second_moment = np.trapezoid(
        combined * radial_grid**3, radial_grid
    ) / np.trapezoid(
        combined * radial_grid, radial_grid
    )
    nlight_scale = 0.5 * 4.0 * sqrt(0.5 * radial_second_moment)
    inner_radius = mode["r0"] / nlight_scale
    inner_sigma = mode["sigma0"] / nlight_scale
    outer_radius = mode["r1"] / nlight_scale
    outer_sigma = mode["sigma1"] / nlight_scale
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
    models.append(
        (
            f"projected\nnLightAFX Index 6 profile + exponential projection ($\\alpha={alpha:.3f}$)",
            *normalize(top, section),
        )
    )

    table_x, table_y, table = read_tabulated(args.tabulated_profile)
    table_centroid, table_d4sigma = tabulated_metrics(table_x, table_y, table)
    table_scale = 0.5 * table_d4sigma
    table_x = table_x / table_scale
    table_y = table_y / table_scale
    row = table[np.argmin(np.abs(table_y - table_centroid[1] / table_scale))]
    lateral_table = np.interp(lateral, table_x, row, left=0.0, right=0.0)
    section = lateral_table[np.newaxis, :] * np.exp(-3.0 * depth[:, np.newaxis] ** 2)
    models.append(
        (
            "projected\ntabulated profile + exponential projection",
            *normalize(table, section),
        )
    )

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
        "(c)  projected\n      superGaussian profile\n      + exponential projection",
        f"(d)  projected\n      nLightAFX Index 6 profile, $\\alpha={alpha:.3f}$",
        "(e)  projected\n      tabulated profile\n      + exponential projection",
    )
    images = []
    panel_axes = []

    for row_index, ((label, top, section), row_label) in enumerate(
        zip(models, row_labels)
    ):
        top_axis = figure.add_subplot(grid[row_index, 0])
        section_axis = figure.add_subplot(grid[row_index, 1])
        top_extent = [-1.5, 1.5, -1.5, 1.5]
        if row_index == 4:
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
        "projected = profile × projection",
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
            xlabel="Normalized depth below beam plane,  $\\zeta/d$",
            ylabel="$p(\\zeta)$",
            xlim=(0.0, 1.5),
            ylim=(0.0, 1.03),
        )
        projection_axis.grid(True, which="both", linestyle="-", alpha=0.1)
        projection_axis.minorticks_on()
        projection_axis.legend(ncol=4, frameon=True, loc="upper right")
        projection_figure.savefig(args.projection_output)
        plt.close(projection_figure)

    plot_profile_metrics(args.metrics_output)


if __name__ == "__main__":
    main()
