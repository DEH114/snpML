"""Generate the JASPAR dataset for the exact 20 QIAGEN test rsIDs.

The script resolves the rsIDs to current GRCh38 SNV alleles, preserves
unsupported/non-SNV inputs in audit files, then reuses jaspar_branch.py to scan
the same allele-specific windows against a local JASPAR PFM collection.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import jaspar_branch
from challenge_batch import (DEFAULT_INPUT, flatten_summary,
                             prepare_qiagen_panel)
from common import Client, utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "jaspar_dataset"
DEFAULT_PFM = ROOT / "data" / "vertebrates.jaspar"
SUMMARY_FIELDS = [
    "rank", "rsid", "resolution_status", "n_snv_alt_alleles", "variant_keys",
    "has_non_snv_alleles", "var_class", "most_severe_consequence",
    "special_note", "protein_status", "n_missense_transcripts",
    "jaspar_status", "n_motif_hits", "reason",
]


def portable_pfm(source: str | Path, out: Path) -> Path:
    source = Path(source)
    if not source.is_file():
        raise FileNotFoundError(
            f"JASPAR PFM file not found: {source}. Download the CORE vertebrate "
            "non-redundant single-batch JASPAR file first."
        )
    destination = out / source.name
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)
    return destination


def run(args: argparse.Namespace) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    if args.flank < 1 or not 0 <= args.threshold <= 1 or args.min_delta < 0:
        raise ValueError("Require flank >= 1, threshold in [0,1], and min-delta >= 0")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be positive")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pfm_path = portable_pfm(args.pfm_file, out)
    client = Client(args.cache, args.offline, args.refresh)
    prepared = prepare_qiagen_panel(
        args.input, out, client, annotate_other=args.annotate_other,
        require_exact=True,
    )

    branch_exit = 1
    if prepared.variants:
        branch_args = argparse.Namespace(
            input=str(prepared.resolved_path), out=str(out), cache=args.cache,
            offline=args.offline, refresh=args.refresh,
            motif_names=None, matrix_ids=None, pfm_file=str(pfm_path),
            flank=args.flank, threshold=args.threshold,
            min_delta=args.min_delta, pseudocount=args.pseudocount,
        )
        branch_exit = jaspar_branch.run(branch_args)

    summary = flatten_summary(prepared.rsids, prepared.statuses,
                              out / "_protein_not_run", out, False, True)
    write_json(out / "dataset_summary.json", summary)
    write_tsv(out / "dataset_summary.tsv", summary, SUMMARY_FIELDS)
    write_json(out / "dataset_manifest.json", {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "dataset_type": "JASPAR allele-specific motif score changes",
        "input": str(args.input),
        "pfm_file": str(pfm_path),
        "n_input_rsids": len(prepared.rsids),
        "n_resolved_snv_alleles": len(prepared.variants),
        "n_rsids_with_resolved_snv": sum(
            status.get("n_snv_alt_alleles", 0) > 0 for status in prepared.statuses
        ),
        "n_unsupported_or_unresolved": sum(
            status["status"] != "resolved_snv" for status in prepared.statuses
        ),
        "branch_exit_status": branch_exit,
        "preparation_source_provenance": list(client.manifest.values()),
        "notes": [
            "Motif scores are heuristic hypotheses, not measured TF binding",
            "No motif hit does not establish absence of regulatory function",
            "The input contains public loci, not an individual's genotype",
        ],
    })
    print(f"JASPAR dataset complete: {len(prepared.rsids)} rsIDs, "
          f"{len(prepared.variants)} resolved SNV alleles. See {out}")
    return branch_exit if prepared.variants else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="Exact 20-rsID QIAGEN TSV in the original order")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument("--cache", default=".cache", help="Shared public-API response cache")
    parser.add_argument("--offline", action="store_true", help="Use cache only")
    parser.add_argument("--refresh", action="store_true", help="Refresh cached requests")
    parser.add_argument("--pfm-file", default=str(DEFAULT_PFM),
                        help="Local JASPAR-format PFM collection")
    parser.add_argument("--flank", type=int, default=40)
    parser.add_argument("--threshold", type=float, default=0.8)
    parser.add_argument("--min-delta", type=float, default=0.1)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--annotate-other", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Keep rsID-level VEP context for unsupported/non-SNV inputs")
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
