"""QIAGEN BioHack 2026 Challenge 1: rsID-first SNV triage + optional analysis branches.

A public rsID is not necessarily one allele or even a single-nucleotide variant.
This runner keeps the supplied input order, expands resolvable SNV ALT alleles,
flags unsupported/ambiguous alleles and preserves an audit trail for each rsID.

Example:
    python challenge_batch.py --out results/qiagen --prepare-only
    python challenge_batch.py --out results/qiagen --vep-only
    python challenge_batch.py --out results/qiagen --pfm-file data/vertebrates.jaspar
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import common
import jaspar_branch
import protein_branch
from common import (Client, APIError, read_inputs, mapping_candidates, sequence,
                    write_json, write_tsv, utc_now)

DEFAULT_INPUT = Path(__file__).resolve().parent / "examples" / "qiagen_20_rsids.tsv"
SPECIAL_NOTES = {
    "rs429358": "APOE epsilon2/epsilon3/epsilon4 interpretation requires rs429358 + rs7412 allele combinations and individual genotype/phase data; this rsID panel has neither.",
    "rs7412": "APOE epsilon2/epsilon3/epsilon4 interpretation requires rs429358 + rs7412 allele combinations and individual genotype/phase data; this rsID panel has neither.",
    "rs113993960": "Known CFTR p.Phe508del is an in-frame deletion; single-base motif/AlphaMissense scoring is not appropriate for the deletion allele.",
    "rs80357713": "Legacy rsID reported as merged into rs80357783 in dbSNP-related resources; verify current identifiers and insertion allele before interpreting.",
}


def classify_mapping(data: dict) -> dict[str, Any]:
    """Describe all main-chromosome GRCh38 alleles; do not infer clinical effects."""
    mappings = [m for m in data.get("mappings", []) if m.get("assembly_name") == "GRCh38"
                and str(m.get("seq_region_name")) in common.CHROMS]
    snv_alleles, non_snv, all_raw = set(), [], []
    for m in mappings:
        alleles = str(m.get("allele_string", "")).upper().split("/")
        all_raw.append({"chrom": m.get("seq_region_name"), "start": m.get("start"),
                        "end": m.get("end"), "allele_string": m.get("allele_string"),
                        "strand": m.get("strand")})
        if not alleles or len(alleles) < 2:
            non_snv.append({"allele_string": m.get("allele_string"), "reason": "no_ALTs"})
            continue
        ref = alleles[0]
        for alt in alleles[1:]:
            if ref in common.DNA and alt in common.DNA and ref != alt and m.get("start") == m.get("end"):
                snv_alleles.add(alt)
            else:
                non_snv.append({"allele_string": m.get("allele_string"),
                                "reference": ref, "alternative": alt,
                                "reason": "non_single_base_allele_or_non_SNV_span"})
    return {"var_class_source": data.get("var_class"),
            "ensembl_most_severe_consequence": data.get("most_severe_consequence"),
            "main_GRCh38_mappings": all_raw,
            "non_snv_alleles": non_snv, "n_candidate_snv_alts": len(snv_alleles)}


def other_vep(client: Client, rsid: str, out: Path) -> dict:
    """Fallback VEP by rsID: consequence information, *not* exact-allele functional scoring."""
    try:
        records = client.ensembl(f"/vep/human/id/{rsid}", {"canonical": 1, "mane": 1})
        if not isinstance(records, list):
            raise ValueError("Expected list of VEP records")
        summaries, transcript_rows = [], []
        for record in records:
            transcripts = record.get("transcript_consequences", []) or []
            genes = sorted({t.get("gene_symbol") or t.get("gene_id") for t in transcripts
                            if t.get("gene_symbol") or t.get("gene_id")})
            summaries.append({"input": record.get("input"),
                              "allele_string": record.get("allele_string"),
                              "most_severe_consequence": record.get("most_severe_consequence"),
                              "genes": genes})
            for t in transcripts:
                transcript_rows.append({"rsid": rsid, "allele_string": record.get("allele_string"),
                    "variant_allele": t.get("variant_allele"), "gene": t.get("gene_symbol"),
                    "transcript": t.get("transcript_id"), "protein_id": t.get("protein_id"),
                    "protein_start": t.get("protein_start"), "protein_end": t.get("protein_end"),
                    "hgvsp": t.get("hgvsp"), "consequences": t.get("consequence_terms", []),
                    "note": "VEP rsID-only record: check exact allele/assembly; no mutant structure predicted"})
        write_json(out / "other_variants_vep" / f"{rsid}.json", records)
        write_tsv(out / "other_variants_vep" / f"{rsid}.transcripts.tsv", transcript_rows,
                  ["rsid", "allele_string", "variant_allele", "gene", "transcript", "protein_id",
                   "protein_start", "protein_end", "hgvsp", "consequences", "note"])
        return {"status": "vep_rsID_annotation_only", "summary": summaries,
                "n_transcript_consequences": len(transcript_rows),
                "caution": "No ALT selected; cannot assert allele-specific phenotype from this summary"}
    except Exception as exc:
        return {"status": "vep_lookup_failed", "error": str(exc)}


def inspect_rsid(client: Client, rsid: str) -> tuple[list, dict]:
    info: dict[str, Any] = {"rsid": rsid, "status": "lookup_failed",
                            "note": "Annotation is conditional on the exact allele",
                            "special_note": SPECIAL_NOTES.get(rsid, "")}
    try:
        data = client.ensembl(f"/variation/human/{rsid}")
    except Exception as exc:
        info["reason"] = str(exc)
        return [], info
    info.update(classify_mapping(data))
    try:
        variants = mapping_candidates(data, rsid)
        for v in variants:
            actual = sequence(client, v.chrom, v.pos, v.pos)
            if actual != v.ref:
                raise ValueError(f"Reference mismatch: {v.key}; reference is {actual}")
        info.update(status="resolved_snv", resolved_variant_keys=[v.key for v in variants],
                    n_snv_alt_alleles=len(variants))
        if info["non_snv_alleles"]:
            info["status"] = "resolved_snv_with_unhandled_alleles"
            info["reason"] = "Other alleles are non-SNV; only single-base alternatives enter SNV branches"
        return variants, info
    except Exception as exc:
        info.update(status="unsupported_or_ambiguous", reason=str(exc), n_snv_alt_alleles=0,
                    resolved_variant_keys=[])
        return [], info


def make_branch_args(args: argparse.Namespace, input_tsv: Path, out: Path) -> argparse.Namespace:
    return argparse.Namespace(input=str(input_tsv), out=str(out), cache=args.cache,
                              offline=args.offline, refresh=args.refresh,
                              all_transcripts=False, no_alphamissense=False,
                              vep_only=args.vep_only,
                              motif_names=None, matrix_ids=None,
                              pfm_file=str(args.pfm_file or (Path(__file__).parent / "examples" / "CEBPA_MA0102.5.jaspar")),
                              flank=40, threshold=args.threshold,
                              min_delta=args.min_delta, pseudocount=0.5)


def flatten_summary(input_ids: list[str], statuses: list[dict], protein_path: Path,
                    jaspar_path: Path, did_protein: bool, did_jaspar: bool) -> list[dict]:
    p_by_rsid = defaultdict(list)
    j_by_rsid = defaultdict(list)
    if did_protein and (protein_path / "variant_status.json").exists():
        import json
        for row in json.loads((protein_path / "variant_status.json").read_text()):
            p_by_rsid[row.get("rsid")].append(row)
    if did_jaspar and (jaspar_path / "variant_status.json").exists():
        import json
        for row in json.loads((jaspar_path / "variant_status.json").read_text()):
            j_by_rsid[row.get("rsid")].append(row)
    rows = []
    for index, (rsid, status) in enumerate(zip(input_ids, statuses), 1):
        ps, js = p_by_rsid.get(rsid, []), j_by_rsid.get(rsid, [])
        protein_stat = "; ".join(sorted({p.get("status", "unknown") for p in ps})) if ps else (
            "not_run" if not did_protein else "branch_failed_or_no_result" if status.get("n_snv_alt_alleles", 0) else "no_supported_snv")
        jaspar_stat = "; ".join(sorted({j.get("status", "unknown") for j in js})) if js else (
            "not_run" if not did_jaspar else "branch_failed_or_no_result" if status.get("n_snv_alt_alleles", 0) else "no_supported_snv")
        rows.append({"rank": index, "rsid": rsid, "resolution_status": status["status"],
                     "n_snv_alt_alleles": status.get("n_snv_alt_alleles", 0),
                     "variant_keys": status.get("resolved_variant_keys", []),
                     "has_non_snv_alleles": bool(status.get("non_snv_alleles")),
                     "var_class": status.get("var_class_source", ""),
                     "most_severe_consequence": status.get("ensembl_most_severe_consequence", ""),
                     "special_note": status.get("special_note", ""),
                     "protein_status": protein_stat,
                     "n_missense_transcripts": sum(p.get("n_missense_transcripts", 0) for p in ps),
                     "jaspar_status": jaspar_stat,
                     "n_motif_hits": sum(j.get("n_hits", 0) for j in js),
                     "reason": status.get("reason", "")})
    return rows


def run(args: argparse.Namespace) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    inputs = read_inputs(args.input)
    rsids = [r.get("rsid", "") for r in inputs]
    if not rsids or any(not common.re.fullmatch(r"rs\d+", x) for x in rsids):
        raise ValueError("Challenge batch accepts rsID input only, one per row")
    if len(set(rsids)) != len(rsids):
        raise ValueError("Duplicate rsIDs found in input; deduplicate the list")
    if args.only_rsids:
        chosen = args.only_rsids.split(",")
        if any(x not in rsids for x in chosen):
            raise ValueError("--only-rsids contains IDs absent from the provided input")
        rsids = [r for r in rsids if r in set(chosen)]
    if args.max_rsids:
        rsids = rsids[:args.max_rsids]
    write_tsv(out / "challenge_input_manifest.tsv", [{"rank": i, "rsid": r}
                                                   for i, r in enumerate(rsids, 1)], ["rank", "rsid"])
    client = Client(args.cache, args.offline, args.refresh)
    all_variants, statuses = [], []
    for i, rsid in enumerate(rsids, 1):
        print(f"[{i}/{len(rsids)}] Resolving {rsid}...", flush=True)
        variants, status = inspect_rsid(client, rsid)
        if args.annotate_other and status["status"] != "resolved_snv":
            status["other_vep"] = other_vep(client, rsid, out)
        statuses.append(status)
        all_variants.extend(variants)
    write_json(out / "resolution_details.json", statuses)
    write_tsv(out / "resolved_snvs.tsv", [v.row() for v in all_variants],
              ["rsid", "chrom", "pos", "ref", "alt", "assembly", "variant_key"])
    write_tsv(out / "unsupported_or_unresolved.tsv",
              [{"rsid": st["rsid"], "status": st["status"], "reason": st.get("reason", ""),
                "non_snv_alleles": st.get("non_snv_alleles", []),
                "vep_lookup": st.get("other_vep", {}).get("status", "not_requested")}
               for st in statuses if st["status"] != "resolved_snv"],
              ["rsid", "status", "reason", "non_snv_alleles", "vep_lookup"])
    did_p = did_j = False
    p_exit = j_exit = None
    branch_errors = []
    branch_input = out / "resolved_snvs.tsv"
    p_dir, j_dir = out / "protein", out / "jaspar"
    if all_variants and not args.prepare_only:
        bargs = make_branch_args(args, branch_input, p_dir)
        if not args.skip_protein:
            print("Running protein branch...", flush=True)
            did_p = True
            try:
                p_exit = protein_branch.run(bargs)
            except Exception as exc:
                branch_errors.append({"branch": "protein", "error": str(exc)})
        if not args.skip_jaspar:
            print("Running JASPAR branch...", flush=True)
            did_j = True
            bargs.out = str(j_dir)
            try:
                j_exit = jaspar_branch.run(bargs)
            except Exception as exc:
                branch_errors.append({"branch": "jaspar", "error": str(exc)})
    write_json(out / "branch_errors.json", branch_errors)
    summary = flatten_summary(rsids, statuses, p_dir, j_dir, did_p, did_j)
    summary_fields = ["rank", "rsid", "resolution_status", "n_snv_alt_alleles", "variant_keys",
                      "has_non_snv_alleles", "var_class", "most_severe_consequence",
                      "special_note", "protein_status", "n_missense_transcripts",
                      "jaspar_status", "n_motif_hits", "reason"]
    write_json(out / "challenge_summary.json", summary)
    write_tsv(out / "challenge_summary.tsv", summary, summary_fields)
    write_json(out / "challenge_manifest.json", {"finished_at_utc": utc_now(),
                "input": str(args.input), "n_input_rsids": len(rsids),
                "n_resolved_snv_alleles": len(all_variants),
                "n_rsids_with_resolved_snv": sum(st["n_snv_alt_alleles"] > 0 for st in statuses),
                "n_rsids_with_unhandled_or_missing_alleles": sum(st["status"] != "resolved_snv" for st in statuses),
                "protein_exit_status": p_exit, "jaspar_exit_status": j_exit,
                "branch_errors": branch_errors,
                "jaspar_motif_scope": str(args.pfm_file) if args.pfm_file else "CEBPA_MA0102.5 ONLY (demonstration)",
                "notes": ["No phenotype is inferred from AlphaFold pLDDT or motif matching alone",
                          "rsID-only inputs expand all currently resolvable single-nucleotide ALT alleles",
                          "indels and ambiguous alleles are not scored by these SNV-only branches"],
                "source_provenance": list(client.manifest.values())})
    print(f"Done. {len(rsids)} rsIDs, {len(all_variants)} resolved SNV alleles. "
          f"See {out / 'challenge_summary.tsv'}")
    return 0 if all_variants else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="A TSV with rsid header or a plain text list of rsIDs")
    parser.add_argument("--out", default="results/qiagen", help="Output directory")
    parser.add_argument("--cache", default=".cache", help="Public API response cache")
    parser.add_argument("--offline", action="store_true", help="Only use previously cached API responses")
    parser.add_argument("--refresh", action="store_true", help="Re-fetch cached requests")
    parser.add_argument("--prepare-only", action="store_true", help="Resolve IDs and stop before analysis")
    parser.add_argument("--vep-only", action="store_true", help="Protein branch: skip UniProt/structure enrichment")
    parser.add_argument("--skip-protein", action="store_true")
    parser.add_argument("--skip-jaspar", action="store_true")
    parser.add_argument("--pfm-file", help="JASPAR-format file; default: single CEBPA example, NOT a full screen")
    parser.add_argument("--threshold", type=float, default=0.8)
    parser.add_argument("--min-delta", type=float, default=0.1)
    parser.add_argument("--max-rsids", type=int, help="Limit to the first N IDs (smoke test)")
    parser.add_argument("--only-rsids", help="Comma-separated subset, e.g. rs334,rs429358")
    parser.add_argument("--annotate-other", action="store_true",
                        help="Also fetch lightweight rsID-level VEP consequences for non-SNV/ambiguous entries")
    args = parser.parse_args()
    if args.offline and args.refresh:
        parser.error("--offline and --refresh are incompatible")
    if args.max_rsids is not None and args.max_rsids < 1:
        parser.error("--max-rsids must be >= 1")
    try:
        return run(args)
    except Exception as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
