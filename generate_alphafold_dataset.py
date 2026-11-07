"""Generate the AlphaFold/protein dataset for the exact 20 QIAGEN test rsIDs.

The script resolves the rsIDs to current GRCh38 SNV alleles, preserves
unsupported/non-SNV inputs in audit files, then reuses protein_branch.py for
VEP, UniProt and existing AlphaFold reference-model annotations.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import protein_branch
from challenge_batch import (DEFAULT_INPUT, flatten_summary,
                             prepare_qiagen_panel)
from common import Client, utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "alpha_fold_dataset"
SUMMARY_FIELDS = [
    "rank", "rsid", "resolution_status", "n_snv_alt_alleles", "variant_keys",
    "has_non_snv_alleles", "var_class", "most_severe_consequence",
    "special_note", "protein_status", "n_missense_transcripts",
    "jaspar_status", "n_motif_hits", "reason",
]


def run(args: argparse.Namespace) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
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
            all_transcripts=args.all_transcripts, vep_only=args.vep_only,
            no_alphamissense=args.no_alphamissense,
        )
        branch_exit = protein_branch.run(branch_args)

    summary = flatten_summary(prepared.rsids, prepared.statuses, out,
                              out / "_jaspar_not_run", True, False)
    write_json(out / "dataset_summary.json", summary)
    write_tsv(out / "dataset_summary.tsv", summary, SUMMARY_FIELDS)
    write_json(out / "dataset_manifest.json", {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "dataset_type": "VEP + UniProt + AlphaFold reference-protein context",
        "input": str(args.input),
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
            "AlphaFold files are existing reference-protein models, not mutant predictions",
            "pLDDT is local model confidence, not variant pathogenicity",
            "The input contains public loci, not an individual's genotype",
        ],
    })
    print(f"AlphaFold dataset complete: {len(prepared.rsids)} rsIDs, "
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
    parser.add_argument("--all-transcripts", action="store_true",
                        help="Deeply enrich every missense transcript")
    parser.add_argument("--vep-only", action="store_true",
                        help="Skip UniProt and AlphaFold downloads")
    parser.add_argument("--no-alphamissense", action="store_true",
                        help="Disable the optional VEP AlphaMissense plugin")
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
