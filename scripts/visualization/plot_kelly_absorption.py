#!/usr/bin/env python3
"""Plot the piecewise Kelly absorption model used by AdditiveFOAM."""

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np

from plot_style import PUBLICATION_STYLE


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def kelly_absorptivity(aspect_ratio, geometry, eta0, eta_min, switch):
    absorptivity = np.full_like(aspect_ratio, eta_min, dtype=float)
    active = aspect_ratio > switch
    a = aspect_ratio[active]
    theta = np.arctan(1.0 / a)
    if geometry == "cone":
        f = (3.0 * np.sin(theta) - np.sin(3.0 * theta)) / 4.0
        g = 1.0 / (1.0 + np.sqrt(1.0 + a**2))
    elif geometry == "cylinder":
        f = (1.0 - np.cos(2.0 * theta)) / 2.0
        g = 1.0 / (2.0 * (1.0 + a))
    else:
        raise ValueError(f"Unknown geometry: {geometry}")
    absorptivity[active] = eta0 * (1.0 + (1.0 - eta0) * (g - f)) / (
        1.0 - (1.0 - eta0) * (1.0 - g)
    )
    return absorptivity


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    eta0 = 0.28
    eta_min = 0.35
    switch = 1.0
    aspect_ratio = np.linspace(0.0, 6.0, 1201)

    plt.rcParams.update(PUBLICATION_STYLE)
    figure, axis = plt.subplots(figsize=(4.5, 3.0), constrained_layout=True)
    axis.axvspan(0.0, switch, color="0.92", zorder=0)
    axis.plot(
        aspect_ratio,
        kelly_absorptivity(aspect_ratio, "cone", eta0, eta_min, switch),
        color="#0072B2",
        label="cone",
    )
    axis.plot(
        aspect_ratio,
        kelly_absorptivity(aspect_ratio, "cylinder", eta0, eta_min, switch),
        color="#D55E00",
        linestyle="--",
        label="cylinder",
    )
    axis.axvline(switch, color="0.25", linestyle=":")
    axis.text(
        0.5 * switch,
        eta_min + 0.025,
        "$\\eta=\\eta_{\\min}$",
        ha="center",
        va="bottom",
        color="0.25",
    )
    axis.annotate(
        "$a_s$",
        xy=(switch, 0.255),
        xytext=(switch + 0.12, 0.285),
        arrowprops={"arrowstyle": "-", "color": "0.25"},
        color="0.25",
    )
    axis.set(
        xlabel="Source aspect ratio,  $a=d_z/\\min(d_x,d_y)$",
        ylabel="Effective absorptivity,  $\\eta(a)$",
        xlim=(0.0, 6.0),
        ylim=(0.25, 1.0),
    )
    axis.set_xticks(np.arange(0.0, 6.1, 1.0))
    axis.set_yticks(np.arange(0.3, 1.01, 0.1))
    axis.grid(True, which="both", linestyle="-", alpha=0.1)
    axis.minorticks_on()
    axis.legend(
        title="$\\eta_0=0.28,\\;\\eta_{\\min}=0.35,\\;a_s=1$",
        frameon=True,
        loc="lower right",
    )
    figure.savefig(args.output)
    plt.close(figure)


if __name__ == "__main__":
    main()
