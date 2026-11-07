"""Run all six rsID dataset generators from one validated panel input.

The child generators run sequentially by default.  They share public-service
response caches and each implements its own host-level throttling, so parallel
processes could both corrupt identical cache writes and exceed an upstream
service's aggregate request rate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from challenge_batch import DEFAULT_INPUT, QIAGEN_20_RSIDS, select_panel_rsids
from common import utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset"
DEFAULT_PFM = ROOT / "data" / "vertebrates.jaspar"
NORMALIZED_INPUT_NAME = "panel_input.tsv"
MANIFEST_NAME = "all_datasets_run_manifest.json"

SOURCE_SPECS = (
    ("clinvar", "generate_clinvar_dataset.py", "clinvar_dataset"),
    ("gwas_catalog", "generate_gwas_catalog_dataset.py", "gwas_catalog_dataset"),
    ("gnomad", "generate_gnomad_dataset.py", "gnomad_dataset"),
    ("gtr", "generate_gtr_dataset.py", "gtr_dataset"),
    ("alphafold", "generate_alphafold_dataset.py", "alpha_fold_dataset"),
    ("jaspar", "generate_jaspar_dataset.py", "jaspar_dataset"),
)


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_dbsnp_resolution(path: str | Path, rsids: list[str]) -> bool:
    """Return true only for an ordered, complete resolution of this panel."""
    path = Path(path)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    if not isinstance(value, list) or len(value) != len(rsids):
        return False
    supplied = []
    for row in value:
        if not isinstance(row, dict):
            return False
        input_rsid = str(row.get("input_rsid", ""))
        current_rsid = str(row.get("current_rsid", ""))
        if not re.fullmatch(r"rs[0-9]+", input_rsid):
            return False
        if not re.fullmatch(r"rs[0-9]+", current_rsid):
            return False
        supplied.append(input_rsid)
    return supplied == rsids and len(set(supplied)) == len(supplied)


def _common_child_arguments(args: argparse.Namespace, normalized_input: Path,
                            source_out: Path, custom_panel: bool) -> list[str]:
    values = [
        "--input", str(normalized_input),
        "--out", str(source_out),
        "--cache", str(Path(args.cache).resolve()),
    ]
    if custom_panel:
        values.append("--allow-custom-panel")
    if args.offline:
        values.append("--offline")
    if args.refresh:
        values.append("--refresh")
    return values


def build_source_command(name: str, script_name: str, source_out: Path,
                         args: argparse.Namespace, normalized_input: Path,
                         custom_panel: bool, clinvar_succeeded: bool = False,
                         rsids: list[str] | None = None) -> list[str]:
    command = [sys.executable, str((ROOT / script_name).resolve())]
    command.extend(_common_child_arguments(
        args, normalized_input, source_out, custom_panel,
    ))
    if name in {"clinvar", "gtr"}:
        command.extend(["--max-records-per-rsid", str(args.max_records_per_rsid)])
    elif name == "gwas_catalog":
        if args.archive:
            command.extend(["--archive", str(Path(args.archive).resolve())])
        if args.no_metadata:
            command.append("--no-metadata")
        command.extend([
            "--max-download-bytes", str(args.max_download_bytes),
            "--max-uncompressed-bytes", str(args.max_uncompressed_bytes),
        ])
        resolution = source_out.parent / "clinvar_dataset" / "dbsnp_resolution.json"
        if (clinvar_succeeded and rsids is not None
                and valid_dbsnp_resolution(resolution, rsids)):
            command.extend(["--dbsnp-resolution", str(resolution.resolve())])
    elif name == "alphafold":
        if args.all_transcripts:
            command.append("--all-transcripts")
        if args.vep_only:
            command.append("--vep-only")
        if args.no_alphamissense:
            command.append("--no-alphamissense")
        command.append("--annotate-other" if args.annotate_other
                       else "--no-annotate-other")
    elif name == "jaspar":
        command.extend([
            "--pfm-file", str(Path(args.pfm_file).resolve()),
            "--flank", str(args.flank),
            "--threshold", str(args.threshold),
            "--min-delta", str(args.min_delta),
            "--pseudocount", str(args.pseudocount),
        ])
        command.append("--annotate-other" if args.annotate_other
                       else "--no-annotate-other")
    return command


def build_combine_command(normalized_input: Path, dataset_root: Path,
                          allow_partial: bool) -> list[str]:
    command = [
        sys.executable,
        str((ROOT / "combine_six_datasets.py").resolve()),
        "--input", str(normalized_input),
        "--dataset-root", str(dataset_root),
        "--out", str(dataset_root / "combined_dataset"),
    ]
    if allow_partial:
        command.append("--allow-partial")
    return command


def _write_log(path: Path, value: str | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value or "", encoding="utf-8")


def execute_step(name: str, command: list[str], output_directory: Path,
                 logs: Path) -> dict[str, Any]:
    stdout_log = logs / f"{name}.stdout.log"
    stderr_log = logs / f"{name}.stderr.log"
    started_at = utc_now()
    started_clock = time.monotonic()
    print(f"Running {name}: {' '.join(command)}", flush=True)
    try:
        result = subprocess.run(
            command,
            cwd=str(ROOT),
            shell=False,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        exit_code: int | None = result.returncode
        stdout, stderr = result.stdout or "", result.stderr or ""
        status = "success" if exit_code == 0 else "failed"
    except OSError as exc:
        exit_code = None
        stdout, stderr = "", f"Unable to start subprocess: {exc}\n"
        status = "failed_to_start"
    _write_log(stdout_log, stdout)
    _write_log(stderr_log, stderr)
    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")
    if stderr:
        print(stderr, end="" if stderr.endswith("\n") else "\n", file=sys.stderr)
    return {
        "name": name,
        "status": status,
        "exit_code": exit_code,
        "started_at_utc": started_at,
        "finished_at_utc": utc_now(),
        "duration_seconds": round(time.monotonic() - started_clock, 6),
        "command": command,
        "output_directory": str(output_directory),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
    }


def skipped_step(name: str, script_name: str, source_out: Path,
                 args: argparse.Namespace, normalized_input: Path,
                 custom_panel: bool, reason: str) -> dict[str, Any]:
    command = build_source_command(
        name, script_name, source_out, args, normalized_input, custom_panel,
    )
    return {
        "name": name,
        "status": "skipped",
        "reason": reason,
        "exit_code": None,
        "started_at_utc": None,
        "finished_at_utc": None,
        "duration_seconds": 0,
        "command": command,
        "output_directory": str(source_out),
        "stdout_log": None,
        "stderr_log": None,
    }


def run(args: argparse.Namespace) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    if args.max_records_per_rsid < 1:
        raise ValueError("--max-records-per-rsid must be positive")
    if args.max_download_bytes < 1 or args.max_uncompressed_bytes < 1:
        raise ValueError("GWAS archive size limits must be positive")

    original_input = Path(args.input).resolve()
    rsids = select_panel_rsids(original_input, require_exact=False)
    custom_panel = tuple(rsids) != QIAGEN_20_RSIDS
    dataset_root = Path(args.out).resolve()
    dataset_root.mkdir(parents=True, exist_ok=True)
    logs = dataset_root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    normalized_input = dataset_root / NORMALIZED_INPUT_NAME
    write_tsv(normalized_input,
              [{"rank": rank, "rsid": rsid} for rank, rsid in enumerate(rsids, 1)],
              ["rank", "rsid"])

    manifest: dict[str, Any] = {
        "started_at_utc": utc_now(),
        "finished_at_utc": None,
        "generator": Path(__file__).name,
        "execution_mode": "sequential_subprocesses",
        "input": str(original_input),
        "input_sha256": file_sha256(original_input),
        "normalized_input": str(normalized_input),
        "normalized_input_sha256": file_sha256(normalized_input),
        "panel_mode": "custom" if custom_panel else "qiagen_20_exact",
        "n_input_rsids": len(rsids),
        "rsids": rsids,
        "dataset_root": str(dataset_root),
        "cache": str(Path(args.cache).resolve()),
        "fail_fast": args.fail_fast,
        "sources": [],
        "combine": None,
        "overall_status": "running",
    }
    manifest_path = dataset_root / MANIFEST_NAME
    write_json(manifest_path, manifest)

    stop = False
    clinvar_succeeded = False
    for name, script_name, folder_name in SOURCE_SPECS:
        source_out = dataset_root / folder_name
        if stop:
            step = skipped_step(
                name, script_name, source_out, args, normalized_input,
                custom_panel, "fail_fast_after_previous_failure",
            )
        else:
            command = build_source_command(
                name, script_name, source_out, args, normalized_input,
                custom_panel, clinvar_succeeded, rsids,
            )
            step = execute_step(name, command, source_out, logs)
            if name == "clinvar":
                resolution = dataset_root / "clinvar_dataset" / "dbsnp_resolution.json"
                clinvar_succeeded = (step["status"] == "success"
                                     and valid_dbsnp_resolution(resolution, rsids))
                step["dbsnp_resolution_valid"] = clinvar_succeeded
            if step["status"] != "success" and args.fail_fast:
                stop = True
        manifest["sources"].append(step)
        write_json(manifest_path, manifest)

    successful_sources = sum(step["status"] == "success"
                             for step in manifest["sources"])
    all_sources_succeeded = successful_sources == len(SOURCE_SPECS)
    if args.skip_combine:
        manifest["combine"] = {
            "name": "combine", "status": "skipped",
            "reason": "skip_combine_requested", "exit_code": None,
        }
    elif all_sources_succeeded or args.allow_partial_combine:
        command = build_combine_command(
            normalized_input, dataset_root,
            allow_partial=not all_sources_succeeded,
        )
        manifest["combine"] = execute_step(
            "combine", command, dataset_root / "combined_dataset", logs,
        )
    else:
        manifest["combine"] = {
            "name": "combine", "status": "skipped",
            "reason": "one_or_more_source_generators_failed",
            "exit_code": None,
        }

    combine_ok = manifest["combine"]["status"] in {"success", "skipped"}
    overall_ok = all_sources_succeeded and combine_ok
    manifest.update({
        "finished_at_utc": utc_now(),
        "n_sources_succeeded": successful_sources,
        "n_sources_failed_or_skipped": len(SOURCE_SPECS) - successful_sources,
        "overall_status": "success" if overall_ok else "failed",
    })
    write_json(manifest_path, manifest)
    print(f"Dataset suite {manifest['overall_status']}: {successful_sources}/"
          f"{len(SOURCE_SPECS)} sources succeeded. See {manifest_path}")
    return 0 if overall_ok else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="TSV with an rsid header or one rsID per line")
    parser.add_argument("--out", default=str(DEFAULT_OUT),
                        help="Root containing the six dataset subfolders")
    parser.add_argument("--cache", default=".cache",
                        help="Shared public-API response cache")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--offline", action="store_true", help="Use cached data only")
    mode.add_argument("--refresh", action="store_true", help="Refresh cached requests")
    parser.add_argument("--fail-fast", action="store_true",
                        help="Stop launching sources after the first failure")
    parser.add_argument("--skip-combine", action="store_true",
                        help="Do not run combine_six_datasets.py")
    parser.add_argument("--allow-partial-combine", action="store_true",
                        help="Run the combiner after source failures using available data")

    parser.add_argument("--max-records-per-rsid", type=int, default=1000,
                        help="ClinVar and GTR per-rsID safety limit")
    parser.add_argument("--archive", help="Local official GWAS Catalog ZIP")
    parser.add_argument("--no-metadata", action="store_true",
                        help="Skip optional GWAS Catalog metadata")
    parser.add_argument("--max-download-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--max-uncompressed-bytes", type=int, default=4_000_000_000)

    parser.add_argument("--all-transcripts", action="store_true",
                        help="AlphaFold branch: enrich every missense transcript")
    parser.add_argument("--vep-only", action="store_true",
                        help="AlphaFold branch: skip UniProt/structure enrichment")
    parser.add_argument("--no-alphamissense", action="store_true",
                        help="AlphaFold branch: disable the VEP AlphaMissense plugin")
    parser.add_argument("--annotate-other", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="AlphaFold/JASPAR: annotate unsupported rsIDs with VEP")

    parser.add_argument("--pfm-file", default=str(DEFAULT_PFM),
                        help="JASPAR-format PFM collection")
    parser.add_argument("--flank", type=int, default=40)
    parser.add_argument("--threshold", type=float, default=0.8)
    parser.add_argument("--min-delta", type=float, default=0.1)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    return parser


def main() -> int:
    parser = build_parser()
    try:
        return run(parser.parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
