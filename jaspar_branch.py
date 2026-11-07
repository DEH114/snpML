"""Compare REF/ALT motif scores at the SAME SNP-overlapping window and strand.
Default: small targeted CEBPA/CEBPB screen, NOT a whole-JASPAR screen.
Scores and gain/loss calls are heuristic, not p-values, binding measurements,
expression changes, target-gene assignments or phenotype predictions.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from urllib.parse import quote, urljoin
import numpy as np
from Bio import motifs as bio_motifs
from common import (Client, add_common_arguments, read_inputs, resolve, revcomp,
                    save_run, sequence, write_json, write_tsv)

JASPAR = "https://jaspar.elixir.no/api/v1"
BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}


@dataclass
class Matrix:
    matrix_id: str
    name: str
    log_odds: np.ndarray
    metadata: dict

    @property
    def length(self) -> int:
        return self.log_odds.shape[1]

    @property
    def minimum(self) -> float:
        return float(self.log_odds.min(axis=0).sum())

    @property
    def maximum(self) -> float:
        return float(self.log_odds.max(axis=0).sum())

    def score(self, seq: str) -> float:
        if len(seq) != self.length or any(b not in BASE_INDEX for b in seq):
            raise ValueError("Motif window has wrong length or non-ACGT bases")
        return float(self.log_odds[[BASE_INDEX[b] for b in seq], np.arange(self.length)].sum())

    def relative(self, score: float) -> float:
        span = self.maximum - self.minimum
        if span <= 1e-12:
            raise ValueError("Uninformative matrix: maximum equals minimum")
        return float((score - self.minimum) / span)


def make_matrix(profile: dict, pseudocount: float = 0.5) -> Matrix:
    """Uniform A/C/G/T background; add pseudocount to each base at each position."""
    if pseudocount <= 0:
        raise ValueError("Pseudocount must be positive")
    pfm = profile.get("pfm", {})
    counts = np.array([pfm[b] for b in BASES], dtype=float)
    if counts.ndim != 2 or counts.shape[0] != 4 or counts.shape[1] < 1:
        raise ValueError("Expected a 4-by-L JASPAR PFM")
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("PFM counts must be finite and nonnegative")
    probabilities = (counts + pseudocount) / (counts.sum(axis=0, keepdims=True) + 4 * pseudocount)
    matrix = Matrix(str(profile["matrix_id"]), str(profile.get("name", "")),
                    np.log2(probabilities / 0.25), profile)
    if matrix.maximum - matrix.minimum <= 1e-12:
        raise ValueError("Uninformative motif")
    return matrix


def search_profiles(client: Client, name: str) -> list[dict]:
    url = f"{JASPAR}/matrix/"
    params = {"search": name, "collection": "CORE", "tax_group": "vertebrates",
              "version": "latest", "page_size": 100, "format": "json"}
    found, visited = [], set()
    while url:
        if url in visited:
            raise ValueError("Repeated pagination URL from JASPAR")
        visited.add(url)
        result = client.json(url, params)
        for profile in result.get("results", []):
            if str(profile.get("name", "")).upper() == name.upper():
                found.append(profile)
        next_url = result.get("next")
        url = urljoin(url, next_url).replace("http://", "https://", 1) if next_url else ""
        params = None
    # Keep the highest returned version per base matrix identifier.
    latest = {}
    for profile in found:
        mid = profile["matrix_id"]
        base, _, ver = mid.partition(".")
        version = int(ver or 0)
        if base not in latest or version > latest[base][0]:
            latest[base] = (version, profile)
    return [x[1] for x in latest.values()]


def load_matrices(client: Client, args: argparse.Namespace) -> tuple[list[Matrix], list[dict]]:
    profiles, errors = [], []
    if args.pfm_file:
        with Path(args.pfm_file).open(encoding="utf-8") as handle:
            for motif in bio_motifs.parse(handle, "jaspar"):
                profiles.append({"matrix_id": motif.matrix_id, "name": motif.name,
                    "pfm": {b: list(map(float, motif.counts[b])) for b in BASES},
                    "source_file": str(args.pfm_file)})
    else:
        ids = list(args.matrix_ids or [])
        names = args.motif_names if args.motif_names is not None else ([] if ids else ["CEBPA", "CEBPB"])
        for name in names:
            try:
                matches = search_profiles(client, name)
                if not matches:
                    raise ValueError(f"No exact TF-name match found for {name}; try --matrix-ids or --pfm-file")
                ids.extend(p["matrix_id"] for p in matches)
            except Exception as exc:
                errors.append({"stage": "motif_lookup", "name": name, "error": str(exc)})
        for mid in dict.fromkeys(ids):
            try:
                profile = client.json(f"{JASPAR}/matrix/{quote(mid, safe='')}/", {"format": "json"})
                profile["source_url"] = f"{JASPAR}/matrix/{mid}/"
                profiles.append(profile)
            except Exception as exc:
                errors.append({"stage": "motif_download", "matrix_id": mid, "error": str(exc)})
    matrices = []
    for profile in profiles:
        try:
            matrices.append(make_matrix(profile, args.pseudocount))
        except Exception as exc:
            errors.append({"stage": "motif_parse", "matrix_id": profile.get("matrix_id"), "error": str(exc)})
    return matrices, errors


def score_change(ref_relative: float, alt_relative: float, threshold: float) -> str:
    if ref_relative < threshold <= alt_relative:
        return "gain_at_heuristic_threshold"
    if alt_relative < threshold <= ref_relative:
        return "loss_at_heuristic_threshold"
    if alt_relative > ref_relative:
        return "stronger_match"
    if alt_relative < ref_relative:
        return "weaker_match"
    return "unchanged_match"


def scan_variant(ref_seq: str, alt_base: str, snp_index: int, region_start: int,
                 matrix: Matrix, threshold: float = 0.8, min_delta: float = 0.1) -> tuple[list[dict], int]:
    """snp_index is 0-based within sequence; all reported genome positions are 1-based inclusive.
    Returns qualifying windows and number of valid strand-specific windows tested.
    """
    ref_seq = ref_seq.upper()
    if not 0 <= snp_index < len(ref_seq) or alt_base not in BASE_INDEX:
        raise ValueError("Invalid SNP index or ALT base")
    if ref_seq[snp_index] == alt_base:
        raise ValueError("ALT equals reference")
    alt_seq = ref_seq[:snp_index] + alt_base + ref_seq[snp_index + 1:]
    length = matrix.length
    start_min, start_max = max(0, snp_index - length + 1), min(snp_index, len(ref_seq) - length)
    results, n_tested = [], 0
    for offset in range(start_min, start_max + 1):
        ref_window, alt_window = ref_seq[offset:offset + length], alt_seq[offset:offset + length]
        if any(b not in BASE_INDEX for b in ref_window + alt_window):
            continue
        for strand in ("+", "-"):
            reference = ref_window if strand == "+" else revcomp(ref_window)
            alternate = alt_window if strand == "+" else revcomp(alt_window)
            ref_score, alt_score = matrix.score(reference), matrix.score(alternate)
            rr, ar = matrix.relative(ref_score), matrix.relative(alt_score)
            n_tested += 1
            if max(rr, ar) < threshold or abs(ar - rr) < min_delta:
                continue
            results.append({"matrix_id": matrix.matrix_id, "tf": matrix.name,
                "motif_start": region_start + offset, "motif_end": region_start + offset + length - 1,
                "strand": strand, "ref_sequence_motif_orientation": reference,
                "alt_sequence_motif_orientation": alternate,
                "snp_position_in_motif_1based": snp_index - offset + 1 if strand == "+" else length - (snp_index - offset),
                "ref_log2_odds": ref_score, "alt_log2_odds": alt_score, "delta_log2_odds": alt_score - ref_score,
                "ref_relative_score": rr, "alt_relative_score": ar, "delta_relative_score": ar - rr,
                "effect": score_change(rr, ar, threshold),
                "target_gene": None, "target_gene_status": "unresolved_requires_independent_evidence",
                "evidence_type": "computational_motif_prediction",
                "score_note": "heuristic scores; not p-values or experimentally established binding effects"})
    return results, n_tested


def run(args: argparse.Namespace) -> int:
    if not 0 <= args.threshold <= 1 or not 0 <= args.min_delta <= 1 or args.flank < 0:
        raise ValueError("threshold and min-delta must be in [0,1]; flank must be nonnegative")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    client = Client(args.cache, args.offline, args.refresh)
    matrices, errors = load_matrices(client, args)
    write_json(out / "motif_profiles.json", [m.metadata for m in matrices])
    if not matrices:
        save_run(out, client, args, errors, {"status": "failed_no_motifs"})
        print("No motifs loaded. See errors.json or use --pfm-file.", file=sys.stderr)
        return 1
    flank = max(args.flank, max(m.length for m in matrices) - 1)
    all_hits, statuses, seen = [], [], set()
    for index, input_row in enumerate(read_inputs(args.input), 1):
        try:
            variants = resolve(client, input_row)
        except Exception as exc:
            errors.append({"input_row": index, "stage": "resolve", "input": input_row, "error": str(exc)})
            continue
        for v in variants:
            if v.key in seen:
                continue
            seen.add(v.key)
            try:
                start, end = max(1, v.pos - flank), v.pos + flank
                ref_seq = sequence(client, v.chrom, start, end)
                snp_index = v.pos - start
                if ref_seq[snp_index] != v.ref:
                    raise ValueError("Reference sequence mismatch")
                context_file = out / "sequences" / f"{v.chrom}_{v.pos}_{v.ref}_{v.alt}.json"
                write_json(context_file, {**v.row(), "start": start, "end": end, "ref_sequence": ref_seq,
                    "alt_sequence": ref_seq[:snp_index] + v.alt + ref_seq[snp_index + 1:]})
                hits, tested = [], 0
                for matrix in matrices:
                    matched, n_tested = scan_variant(ref_seq, v.alt, snp_index, start, matrix,
                                                     args.threshold, args.min_delta)
                    hits.extend({**v.row(), **h} for h in matched)
                    tested += n_tested
                hits.sort(key=lambda h: (-abs(h["delta_relative_score"]), h["matrix_id"], h["motif_start"], h["strand"]))
                all_hits.extend(hits)
                statuses.append({**v.row(), "status": "hits" if hits else "no_hits_under_heuristic",
                                 "n_motifs": len(matrices), "n_scored_windows": tested, "n_hits": len(hits),
                                 "note": "No hit does not establish absence of regulatory effect"})
            except Exception as exc:
                statuses.append({**v.row(), "status": "analysis_error", "error": str(exc)})
                errors.append({"variant_key": v.key, "stage": "motif_analysis", "error": str(exc)})
    write_json(out / "motif_results.json", all_hits)
    write_json(out / "variant_status.json", statuses)
    fields = ["variant_key", "rsid", "chrom", "pos", "ref", "alt", "matrix_id", "tf", "motif_start", "motif_end", "strand",
              "ref_sequence_motif_orientation", "alt_sequence_motif_orientation", "snp_position_in_motif_1based",
              "ref_log2_odds", "alt_log2_odds", "delta_log2_odds", "ref_relative_score", "alt_relative_score",
              "delta_relative_score", "effect", "target_gene", "target_gene_status", "evidence_type"]
    write_tsv(out / "motif_results.tsv", all_hits, fields)
    save_run(out, client, args, errors, {"n_variants": len(seen), "n_motifs": len(matrices), "n_hits": len(all_hits),
        "n_errors": len(errors), "effective_flank": flank, "pseudocount_per_base": args.pseudocount,
        "background": {b: 0.25 for b in BASES}, "multiple_testing": "not performed; no inferential p-values claimed",
        "screen_scope": "exactly the matrix IDs listed in motif_profiles.json; not necessarily comprehensive"})
    print(f"Saved {len(all_hits)} qualifying windows from {len(matrices)} motifs; {len(errors)} errors to {out}")
    return 0 if any(s.get("status") in {"hits", "no_hits_under_heuristic"} for s in statuses) else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument("--motif-names", nargs="+", help="Exact TF names; default CEBPA CEBPB")
    parser.add_argument("--matrix-ids", nargs="+", help="Versioned JASPAR matrix identifiers")
    parser.add_argument("--pfm-file", help="Local JASPAR-format PFM file; overrides online motif selection")
    parser.add_argument("--flank", type=int, default=40)
    parser.add_argument("--threshold", type=float, default=0.8, help="Heuristic minimum REF or ALT relative score")
    parser.add_argument("--min-delta", type=float, default=0.1, help="Heuristic minimum absolute relative-score change")
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()
    try:
        return run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
