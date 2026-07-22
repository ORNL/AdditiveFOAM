#!/usr/bin/env python3
"""Plot the response curves and final fit from a calibration campaign."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

from plot_style import PUBLICATION_STYLE


COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00")
NORMALIZED_DEPTH_COLUMN = "x_depth_over_half_d4sigma"
EXPECTED_X_DEFINITION = "d / (D4sigma / 2)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot AdditiveFOAM projected-source calibration outputs."
    )
    parser.add_argument("--experiments", type=Path, required=True)
    parser.add_argument("--simulations", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--fit", type=Path, required=True)
    parser.add_argument("--responses-output", type=Path, required=True)
    parser.add_argument("--fit-output", type=Path, required=True)
    return parser.parse_args()


def load_yaml(path: Path):
    with path.open(encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def load_summary(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


def validate_normalization(summary, fit) -> None:
    if not summary or NORMALIZED_DEPTH_COLUMN not in summary[0]:
        raise ValueError(
            "Calibration summary must contain "
            f"{NORMALIZED_DEPTH_COLUMN!r} for the 2sigma normalization."
        )
    if fit.get("x_definition") != EXPECTED_X_DEFINITION:
        raise ValueError(
            "Calibration fit does not declare the expected normalization "
            f"{EXPECTED_X_DEFINITION!r}."
        )


def condition_key(parameters: dict) -> tuple:
    names = ("Power_W", "Speed_mm_s", "Spot_Size_microns")
    return tuple((name, float(parameters[name])) for name in names)


def response_curve(trial_values, depths):
    trial_values = np.asarray(trial_values, dtype=float)
    depths = np.asarray(depths, dtype=float)
    order = np.argsort(trial_values)
    trial_values = trial_values[order]
    depths = depths[order]
    dense_trials = np.linspace(trial_values.min(), trial_values.max(), 500)
    if len(np.unique(trial_values)) >= 3:
        dense_depths = PchipInterpolator(trial_values, depths)(dense_trials)
    else:
        dense_depths = np.interp(dense_trials, trial_values, depths)
    return trial_values, depths, dense_trials, dense_depths


def plot_responses(experiments, simulations, summary, output: Path) -> None:
    summary_by_condition = {condition_key(row): row for row in summary}
    experiments_by_condition = {
        condition_key(row["parameters"]): row for row in experiments
    }
    figure, axes = plt.subplots(
        3,
        2,
        figsize=(6.5, 6.2),
        constrained_layout=True,
        sharex=True,
    )
    axes = axes.ravel()

    for index, simulation in enumerate(simulations):
        axis = axes[index]
        parameters = simulation["parameters"]
        power = float(parameters["Power_W"])
        key = condition_key(parameters)
        summary_row = summary_by_condition[key]
        measured_depths = np.asarray(
            experiments_by_condition[key]["Measured_Depth_microns"], dtype=float
        )
        trials, depths, dense_trials, dense_depths = response_curve(
            simulation["n"], simulation["Simulated_Depth_microns"]
        )

        measured_mean = summary_row["mean_depth_microns"]
        spot_size = summary_row["Spot_Size_microns"]
        calibrated_n = summary_row["calibrated_n"]
        calibrated_depth = np.interp(calibrated_n, dense_trials, dense_depths)

        axis.plot(dense_trials, dense_depths, color=COLORS[index], zorder=2)
        axis.scatter(
            trials,
            depths,
            s=16,
            facecolor="white",
            edgecolor="0.15",
            linewidth=0.7,
            zorder=3,
        )
        axis.axhline(
            measured_mean,
            color="#6A3D9A",
            linestyle="--",
            linewidth=1.0,
            zorder=1,
        )
        axis.axhspan(
            measured_depths.min(),
            measured_depths.max(),
            color="#6A3D9A",
            alpha=0.10,
            zorder=0,
        )
        axis.scatter(
            [calibrated_n],
            [calibrated_depth],
            marker="*",
            s=48,
            color="#007A53",
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        axis.set_title(
            f"{power:g} W   ($D_{{4\\sigma}}={spot_size:g}$ µm)", pad=4
        )
        axis.set_xlim(float(trials.min()), float(trials.max()))
        axis.set_ylim(bottom=0.0)
        axis.set_xticks([0, 3, 6, 9])
        axis.grid(True, linestyle="-", alpha=0.1)
        axis.minorticks_on()

    for axis in axes[::2]:
        axis.set_ylabel("Maximum liquidus depth (µm)")
    for axis in axes[4:]:
        axis.set_xlabel("Local shape parameter,  $n=B$")

    legend_axis = axes[-1]
    legend_axis.axis("off")
    legend_axis.legend(
        handles=(
            Line2D([], [], color=COLORS[0], label="PCHIP response"),
            Line2D(
                [],
                [],
                marker="o",
                markerfacecolor="white",
                markeredgecolor="0.15",
                linestyle="none",
                label="AdditiveFOAM trial",
            ),
            Line2D(
                [],
                [],
                color="#6A3D9A",
                linestyle="--",
                label="Mean measured depth",
            ),
            Patch(
                facecolor="#6A3D9A",
                alpha=0.10,
                label="Measured replicate range",
            ),
            Line2D(
                [],
                [],
                marker="*",
                color="#007A53",
                markeredgecolor="white",
                linestyle="none",
                markersize=9,
                label="Local posterior mode",
            ),
        ),
        loc="center",
        frameon=True,
    )
    legend_axis.text(
        0.5,
        0.16,
        "Scan speed: 500 mm/s",
        ha="center",
        va="center",
        transform=legend_axis.transAxes,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def weighted_initial_guess(z_obs, n_obs, sigma):
    design = np.column_stack([z_obs, np.ones_like(z_obs)])
    weighted_design = design / sigma[:, None]
    weighted_values = n_obs / sigma
    theta, *_ = np.linalg.lstsq(weighted_design, weighted_values, rcond=None)
    return theta


def robust_fit(z_obs, n_obs, sigma, loss, f_scale):
    initial = weighted_initial_guess(z_obs, n_obs, sigma)

    def residual(theta):
        return (n_obs - (theta[0] * z_obs + theta[1])) / sigma

    return least_squares(
        residual,
        x0=initial,
        loss=loss,
        f_scale=f_scale,
    ).x


def bootstrap_band(summary, states, fit, x_range):
    state_by_key = {
        condition_key(state["parameters"]): state for state in states
    }
    z_obs = np.log2([row[NORMALIZED_DEPTH_COLUMN] for row in summary])
    n_mode = np.asarray([row["calibrated_n"] for row in summary], dtype=float)
    sigma = np.asarray([row["calibrated_n_std"] for row in summary], dtype=float)
    sigma = np.maximum(sigma, float(fit.get("local_n_sigma_floor", 1.0e-6)))
    posterior_samples = []
    for row in summary:
        parameters = {
            "Power_W": row["Power_W"],
            "Speed_mm_s": row["Speed_mm_s"],
            "Spot_Size_microns": row["Spot_Size_microns"],
        }
        state = state_by_key[condition_key(parameters)]
        posterior_samples.append(
            np.asarray(state["posterior_n_samples"], dtype=float)
        )

    rng = np.random.default_rng(int(fit.get("bootstrap_random_seed", 12345)))
    requested = int(fit.get("bootstrap_samples_requested", 1000))
    loss = str(fit.get("loss", "soft_l1"))
    f_scale = float(fit.get("f_scale", 1.0))
    z_range = np.log2(x_range)
    curves = []
    for _ in range(requested):
        sampled_n = np.asarray(
            [rng.choice(samples) for samples in posterior_samples], dtype=float
        )
        sampled_n = np.clip(sampled_n, 0.0, 9.0)
        theta = robust_fit(z_obs, sampled_n, sigma, loss, f_scale)
        curves.append(np.clip(theta[0] * z_range + theta[1], 0.0, 9.0))
    return np.quantile(np.asarray(curves), [0.025, 0.5, 0.975], axis=0)


def plot_fit(summary, states, fit, output: Path) -> None:
    x_obs = np.asarray(
        [row[NORMALIZED_DEPTH_COLUMN] for row in summary], dtype=float
    )
    n_obs = np.asarray([row["calibrated_n"] for row in summary], dtype=float)
    lower = np.asarray([row["calibrated_n_95_low"] for row in summary])
    upper = np.asarray([row["calibrated_n_95_high"] for row in summary])
    powers = np.asarray([row["Power_W"] for row in summary])
    x_range = np.geomspace(0.95 * x_obs.min(), 1.05 * x_obs.max(), 300)
    band_lower, band_median, band_upper = bootstrap_band(
        summary, states, fit, x_range
    )
    fitted = np.clip(
        float(fit["A"]) * np.log2(x_range) + float(fit["B"]), 0.0, 9.0
    )

    figure, axis = plt.subplots(figsize=(4.5, 3.25), constrained_layout=True)
    axis.fill_between(
        x_range,
        band_lower,
        band_upper,
        color="#E69F00",
        alpha=0.22,
        label="Empirical 95% fit interval",
        zorder=1,
    )
    axis.plot(
        x_range,
        fitted,
        color="#007A53",
        label=(
            "$n=A\\log_2(x)+B$\n"
            f"$A={float(fit['A']):.3f}$, $B={float(fit['B']):.3f}$"
        ),
        zorder=2,
    )
    axis.plot(
        x_range,
        band_median,
        color="#E69F00",
        linestyle=":",
        linewidth=0.9,
        label="Bootstrap median",
        zorder=2,
    )
    for index, (x_value, n_value, low, high, power) in enumerate(
        zip(x_obs, n_obs, lower, upper, powers)
    ):
        axis.errorbar(
            x_value,
            n_value,
            yerr=[[n_value - low], [high - n_value]],
            fmt="o",
            color=COLORS[index],
            markeredgecolor="white",
            markeredgewidth=0.6,
            capsize=2.5,
            zorder=4,
        )
        axis.annotate(
            f"{power:g} W",
            (x_value, n_value),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    axis.set_xscale("log", base=2)
    axis.set(
        xlabel="Normalized measured depth,  $x=d/(2\\sigma)$",
        ylabel="Local shape parameter,  $n$",
        xlim=(x_range.min(), x_range.max()),
        ylim=(0.0, 7.0),
    )
    axis.grid(True, which="both", linestyle="-", alpha=0.1)
    axis.minorticks_on()
    axis.legend(frameon=True, loc="upper left")

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    plt.rcParams.update(PUBLICATION_STYLE)
    simulations = load_yaml(args.simulations)
    experiments = load_yaml(args.experiments)
    summary = load_summary(args.summary)
    states = load_yaml(args.state)
    fit = load_yaml(args.fit)
    validate_normalization(summary, fit)
    plot_responses(experiments, simulations, summary, args.responses_output)
    plot_fit(summary, states, fit, args.fit_output)


if __name__ == "__main__":
    main()
