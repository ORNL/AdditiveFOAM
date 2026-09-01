#!/usr/bin/env python3
"""Regenerate every simulation-derived documentation visualization.

The OpenFOAM cases are inputs to this driver and remain outside the
documentation repository. Each case must already be run and reconstructed.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_OUTPUT_DIR = REPOSITORY_ROOT / "assets/images/visualizations"

CASE_DIRECTORY_NAMES = {
    "baseline": "AMB2018-02-B",
    "multi_beam": "multiBeam",
    "nlight": "nLightAFX",
    "tabulated": "tabulated",
    "amr": "AMB2018-02-B-AMR",
    "multi_layer": "multiLayerPBF",
}

RENDER_GROUPS = (
    "amb2018",
    "multi-beam",
    "nlight-afx",
    "tabulated",
    "multi-layer",
    "amr",
    "multi-layer-fields",
)

EXPECTED_ASSETS = (
    "amb2018-cet.png",
    "amb2018-melt-pool-dimensions.png",
    "amb2018-power.png",
    "amb2018-temperature.mp4",
    "amr-refinement.mp4",
    "amr-refinement.png",
    "heat-source-calibration-fit.png",
    "heat-source-calibration-responses.png",
    "heat-source-models.png",
    "heat-source-projection.png",
    "kelly-absorption.png",
    "multi-beam-temperature.mp4",
    "multi-beam-temperature.png",
    "multi-layer-fields.png",
    "multi-layer-temperature.mp4",
    "multi-layer-temperature.png",
    "nlight-afx-temperature.mp4",
    "nlight-afx-temperature.png",
    "profile-metrics.png",
    "quick-start-temperature.png",
    "tabulated-temperature.mp4",
    "tabulated-temperature.png",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate all AdditiveFOAM documentation figures and videos "
            "from completed tutorial cases and calibration outputs."
        )
    )
    parser.add_argument(
        "--cases-root",
        type=Path,
        default=os.environ.get("ADDITIVEFOAM_DOC_CASES"),
        help=(
            "Directory containing the standard case subdirectories. Defaults "
            "to ADDITIVEFOAM_DOC_CASES."
        ),
    )
    parser.add_argument("--baseline-case", type=Path)
    parser.add_argument(
        "--amb2018-data-case",
        type=Path,
        help=(
            "AMB2018 case containing log.additiveFoam, melt-pool CSVs, and "
            "solidification data for the three quantitative plots. Defaults "
            "to --baseline-case."
        ),
    )
    parser.add_argument("--multi-beam-case", type=Path)
    parser.add_argument("--nlight-case", type=Path)
    parser.add_argument("--tabulated-case", type=Path)
    parser.add_argument(
        "--nlight-config",
        type=Path,
        default=(
            Path(os.environ["ADDITIVEFOAM_ETC"])
            / "heatSources"
            / "nLightAFX-1000.cfg"
            if os.environ.get("ADDITIVEFOAM_ETC")
            else None
        ),
        help=(
            "Shared nLight AFX mode configuration. Defaults to "
            "$ADDITIVEFOAM_ETC/heatSources/nLightAFX-1000.cfg."
        ),
    )
    parser.add_argument("--amr-case", type=Path)
    parser.add_argument("--multi-layer-case", type=Path)
    parser.add_argument(
        "--calibration-campaign",
        type=Path,
        default=os.environ.get("ADDITIVEFOAM_CALIBRATION_CAMPAIGN"),
        help=(
            "Completed heat-source calibration campaign. Defaults to "
            "ADDITIVEFOAM_CALIBRATION_CAMPAIGN or "
            "<cases-root>/heatSourceCalibration/campaign."
        ),
    )
    parser.add_argument(
        "--calibration-experiments",
        type=Path,
        default=os.environ.get("ADDITIVEFOAM_CALIBRATION_EXPERIMENTS"),
        help=(
            "Experimental YAML used by the calibration campaign. Defaults to "
            "ADDITIVEFOAM_CALIBRATION_EXPERIMENTS or experiments.yml beside "
            "the campaign directory."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--pvbatch",
        default="pvbatch",
        help="ParaView batch executable (default: pvbatch on PATH).",
    )
    parser.add_argument(
        "--paraview-pythonpath",
        type=Path,
        default=Path("/usr/lib/python3/dist-packages"),
        help="Directory containing the system ParaView Python modules.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands without executing them.",
    )
    return parser.parse_args()


def command_path(name: str, dry_run: bool) -> str:
    path = shutil.which(name)
    if path:
        return path
    if dry_run:
        return name
    raise SystemExit(
        "Required command {!r} was not found on PATH. Source the OpenFOAM "
        "and AdditiveFOAM environments before running this script.".format(name)
    )


def validate_case(name: str, path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise SystemExit("{} case directory does not exist: {}".format(name, resolved))
    marker = resolved / "case.foam"
    if not marker.is_file():
        raise SystemExit(
            "{} is missing. Reconstruct the case and create this empty ParaView "
            "marker before rendering.".format(marker)
        )
    return resolved


def validate_calibration_campaign(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    required = (
        resolved / "simulations.yml",
        resolved / "calibration_state.yml",
        resolved / "calibration_fit.yml",
        resolved / "reports" / "calibration_summary.csv",
    )
    missing = [str(item) for item in required if not item.is_file()]
    if missing:
        raise SystemExit(
            "Calibration campaign is incomplete; missing: {}".format(
                ", ".join(missing)
            )
        )
    return resolved


def validate_amb2018_data(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    required = (
        resolved / "log.additiveFoam",
        resolved / "solidificationData" / "solidification-data.csv",
    )
    missing = [str(item) for item in required if not item.is_file()]
    dimensions = list(
        (resolved / "postProcessing" / "meltPoolDimensions").glob("*.csv")
    )
    if not dimensions:
        missing.append(
            str(resolved / "postProcessing" / "meltPoolDimensions" / "*.csv")
        )
    if missing:
        raise SystemExit(
            "AMB2018 quantitative-plot data are incomplete; missing: {}. "
            "Run the documented function objects or pass --amb2018-data-case."
            .format(", ".join(missing))
        )
    return resolved


def select_case(
    label: str,
    override: Path | None,
    cases_root: Path | None,
    directory_name: str,
) -> Path:
    if override is not None:
        return validate_case(label, override)
    if cases_root is None:
        raise SystemExit(
            "Set ADDITIVEFOAM_DOC_CASES, pass --cases-root, or provide the "
            "individual --*-case arguments."
        )
    return validate_case(label, cases_root / directory_name)


def run(command: list[str], environment: dict[str, str], dry_run: bool) -> None:
    print("+ {}".format(" ".join(command)), flush=True)
    if not dry_run:
        subprocess.run(command, check=True, env=environment)


def main() -> None:
    args = parse_args()
    cases_root = (
        args.cases_root.expanduser().resolve() if args.cases_root is not None else None
    )
    cases = {
        "baseline": select_case(
            "Baseline",
            args.baseline_case,
            cases_root,
            CASE_DIRECTORY_NAMES["baseline"],
        ),
        "multi_beam": select_case(
            "Multi-beam",
            args.multi_beam_case,
            cases_root,
            CASE_DIRECTORY_NAMES["multi_beam"],
        ),
        "nlight": select_case(
            "nLight AFX",
            args.nlight_case,
            cases_root,
            CASE_DIRECTORY_NAMES["nlight"],
        ),
        "tabulated": select_case(
            "Tabulated",
            args.tabulated_case,
            cases_root,
            CASE_DIRECTORY_NAMES["tabulated"],
        ),
        "amr": select_case(
            "AMR", args.amr_case, cases_root, CASE_DIRECTORY_NAMES["amr"]
        ),
        "multi_layer": select_case(
            "Multi-layer",
            args.multi_layer_case,
            cases_root,
            CASE_DIRECTORY_NAMES["multi_layer"],
        ),
    }
    amb2018_data = validate_amb2018_data(
        args.amb2018_data_case
        if args.amb2018_data_case is not None
        else cases["baseline"]
    )
    calibration_campaign_input = args.calibration_campaign
    if calibration_campaign_input is None:
        if cases_root is None:
            raise SystemExit(
                "Set ADDITIVEFOAM_CALIBRATION_CAMPAIGN, pass "
                "--calibration-campaign, or use --cases-root with a completed "
                "heatSourceCalibration/campaign directory."
            )
        calibration_campaign_input = (
            cases_root / "heatSourceCalibration" / "campaign"
        )
    calibration_campaign = validate_calibration_campaign(
        calibration_campaign_input
    )
    calibration_experiments = args.calibration_experiments
    if calibration_experiments is None:
        calibration_experiments = calibration_campaign.parent / "experiments.yml"
    calibration_experiments = calibration_experiments.expanduser().resolve()
    if not calibration_experiments.is_file():
        raise SystemExit(
            "Calibration experiments file does not exist: {}".format(
                calibration_experiments
            )
        )
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.nlight_config is None:
        raise SystemExit(
            "Source the AdditiveFOAM environment or pass --nlight-config."
        )
    nlight_config = args.nlight_config.expanduser().resolve()
    if not nlight_config.is_file():
        raise SystemExit(
            "nLight AFX configuration does not exist: {}".format(nlight_config)
        )

    environment = os.environ.copy()
    environment["QT_QPA_PLATFORM"] = "offscreen"
    # OpenFOAM can prepend ParaView's bundled Python packages globally. Keep
    # those packages out of Matplotlib/pandas subprocesses and add them back
    # only for pvbatch below.
    environment.pop("PYTHONPATH", None)
    render_environment = environment.copy()
    pythonpath = args.paraview_pythonpath.expanduser().resolve()
    if pythonpath.is_dir():
        existing = render_environment.get("PYTHONPATH")
        render_environment["PYTHONPATH"] = (
            str(pythonpath) if not existing else "{}:{}".format(pythonpath, existing)
        )

    pbatch = command_path(args.pvbatch, args.dry_run)
    ffmpeg = command_path("ffmpeg", args.dry_run)
    del ffmpeg  # render_documentation.py invokes it through PATH.
    plot_power = command_path("plotPower", args.dry_run)
    plot_dimensions = command_path("plotDimensions", args.dry_run)
    plot_cet = command_path("plotCET", args.dry_run)

    quantitative_plots = (
        (plot_power, "amb2018-power.png"),
        (plot_dimensions, "amb2018-melt-pool-dimensions.png"),
        (plot_cet, "amb2018-cet.png"),
    )
    for executable, filename in quantitative_plots:
        run(
            [executable, str(amb2018_data), "-o", str(output_dir / filename)],
            environment,
            args.dry_run,
        )

    render_base = [
        pbatch,
        "--force-offscreen-rendering",
        str(SCRIPT_DIR / "render_documentation.py"),
        "--baseline-case",
        str(cases["baseline"]),
        "--multi-beam-case",
        str(cases["multi_beam"]),
        "--nlight-case",
        str(cases["nlight"]),
        "--tabulated-case",
        str(cases["tabulated"]),
        "--amr-case",
        str(cases["amr"]),
        "--multi-layer-case",
        str(cases["multi_layer"]),
        "--output-dir",
        str(output_dir),
    ]
    for group in RENDER_GROUPS:
        run(render_base + ["--only", group], render_environment, args.dry_run)

    run(
        [
            sys.executable,
            str(SCRIPT_DIR / "plot_heat_source_models.py"),
            "--tabulated-profile",
            str(cases["tabulated"] / "constant/beamProfile.txt"),
            "--nlight-config",
            str(nlight_config),
            "--output",
            str(output_dir / "heat-source-models.png"),
            "--projection-output",
            str(output_dir / "heat-source-projection.png"),
            "--metrics-output",
            str(output_dir / "profile-metrics.png"),
        ],
        environment,
        args.dry_run,
    )
    run(
        [
            sys.executable,
            str(SCRIPT_DIR / "plot_kelly_absorption.py"),
            "--output",
            str(output_dir / "kelly-absorption.png"),
        ],
        environment,
        args.dry_run,
    )
    run(
        [
            sys.executable,
            str(SCRIPT_DIR / "plot_heat_source_calibration.py"),
            "--experiments",
            str(calibration_experiments),
            "--simulations",
            str(calibration_campaign / "simulations.yml"),
            "--summary",
            str(calibration_campaign / "reports" / "calibration_summary.csv"),
            "--state",
            str(calibration_campaign / "calibration_state.yml"),
            "--fit",
            str(calibration_campaign / "calibration_fit.yml"),
            "--responses-output",
            str(output_dir / "heat-source-calibration-responses.png"),
            "--fit-output",
            str(output_dir / "heat-source-calibration-fit.png"),
        ],
        environment,
        args.dry_run,
    )

    if args.dry_run:
        return
    missing = [name for name in EXPECTED_ASSETS if not (output_dir / name).is_file()]
    if missing:
        raise SystemExit("Regeneration completed with missing assets: {}".format(", ".join(missing)))
    print(
        "Regenerated and verified {} documentation assets in {}".format(
            len(EXPECTED_ASSETS), output_dir
        )
    )


if __name__ == "__main__":
    main()
