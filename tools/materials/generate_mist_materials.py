#!/usr/bin/env python3
"""Generate AdditiveFOAM material inputs from Myna Mist material data."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MYNA_URL = "https://github.com/ORNL-MDF/Myna.git"
DEFAULT_MIST_URL = "https://github.com/ORNL-MDF/mist.git"


def run(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, cwd=cwd, text=True).strip()


def load_manifest(output_dir: Path) -> dict:
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    with manifest_path.open(encoding="utf-8") as f:
        return json.load(f)


def clone_repo(url: str, ref: str, destination: Path) -> str:
    run(["git", "clone", "--quiet", url, str(destination)])
    run(["git", "checkout", "--quiet", ref], cwd=destination)
    return run(["git", "rev-parse", "HEAD"], cwd=destination)


def install_mist(mist_url: str, mist_ref: str) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        mist_dir = Path(tmp) / "mist"
        resolved_ref = clone_repo(mist_url, mist_ref, mist_dir)
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--quiet", str(mist_dir)]
        )
        return resolved_ref


def import_mistlib(install: bool, mist_url: str, mist_ref: str) -> tuple[object, str]:
    resolved_ref = mist_ref
    if install:
        resolved_ref = install_mist(mist_url, mist_ref)

    try:
        import mistlib as mist  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "mistlib is not importable. Re-run with --install-mist or install Mist."
        ) from exc

    return mist, resolved_ref


def material_files(myna_dir: Path) -> list[Path]:
    data_dir = myna_dir / "src" / "myna" / "mist_material_data"
    if not data_dir.is_dir():
        raise SystemExit(f"Myna material data directory not found: {data_dir}")
    return sorted(data_dir.glob("*.json"))


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def generate_material(
    mist: object,
    material_json: Path,
    myna_dir: Path,
    output_dir: Path,
) -> dict:
    material_name = material_json.stem
    material_dir = output_dir / material_name
    material_dir.mkdir(parents=True, exist_ok=True)

    mat = mist.core.MaterialInformation(str(material_json))
    transport_file = material_dir / "transportProperties"
    thermo_file = material_dir / "thermoPath"
    mat.write_additivefoam_input(
        transport_file=str(transport_file),
        thermo_file=str(thermo_file),
    )

    return {
        "name": material_name,
        "source": str(material_json.relative_to(myna_dir)),
        "files": [
            display_path(transport_file),
            display_path(thermo_file),
        ],
    }


def write_manifest(
    output_dir: Path,
    myna_url: str,
    myna_ref: str,
    myna_sha: str,
    mist_url: str,
    mist_ref: str,
    mist_sha: str,
    materials: list[dict],
) -> None:
    manifest = {
        "generated_by": "tools/materials/generate_mist_materials.py",
        "sources": {
            "myna": {
                "url": myna_url,
                "ref": myna_ref,
                "sha": myna_sha,
                "material_data_path": "src/myna/mist_material_data",
            },
            "mist": {
                "url": mist_url,
                "ref": mist_ref,
                "sha": mist_sha,
                "package": "mistlib",
            },
        },
        "materials": materials,
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "materials"))
    parser.add_argument("--myna-url", default=DEFAULT_MYNA_URL)
    parser.add_argument("--mist-url", default=DEFAULT_MIST_URL)
    parser.add_argument("--myna-ref", default="main")
    parser.add_argument("--mist-ref", default="main")
    parser.add_argument("--myna-dir", default=None)
    parser.add_argument(
        "--from-manifest",
        action="store_true",
        help="Use source URLs and SHAs from materials/manifest.json.",
    )
    parser.add_argument(
        "--install-mist",
        action="store_true",
        help="Install Mist from the selected Git ref before generation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    myna_checkout_ref = args.myna_ref
    mist_install_ref = args.mist_ref

    if args.from_manifest:
        manifest = load_manifest(output_dir)
        sources = manifest.get("sources", {})
        myna = sources.get("myna", {})
        mist = sources.get("mist", {})
        args.myna_url = myna.get("url", args.myna_url)
        args.myna_ref = myna.get("ref", args.myna_ref)
        myna_checkout_ref = myna.get("sha", args.myna_ref)
        args.mist_url = mist.get("url", args.mist_url)
        args.mist_ref = mist.get("ref", args.mist_ref)
        mist_install_ref = mist.get("sha", args.mist_ref)

    mist_module, mist_sha = import_mistlib(
        args.install_mist,
        args.mist_url,
        mist_install_ref,
    )

    with tempfile.TemporaryDirectory() as tmp:
        if args.myna_dir:
            myna_dir = Path(args.myna_dir).resolve()
            myna_sha = run(["git", "rev-parse", "HEAD"], cwd=myna_dir)
        else:
            myna_dir = Path(tmp) / "Myna"
            myna_sha = clone_repo(args.myna_url, myna_checkout_ref, myna_dir)

        for path in output_dir.iterdir():
            if path.name == "manifest.json":
                continue
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()

        materials = [
            generate_material(mist_module, material_json, myna_dir, output_dir)
            for material_json in material_files(myna_dir)
        ]

    write_manifest(
        output_dir,
        args.myna_url,
        args.myna_ref,
        myna_sha,
        args.mist_url,
        args.mist_ref,
        mist_sha,
        materials,
    )


if __name__ == "__main__":
    main()
