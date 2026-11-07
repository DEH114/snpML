"""Exact-allele VEP -> UniProt -> existing AlphaFold model context.
This script DOES NOT predict mutant structures, stability or a patient's phenotype.
Run: python protein_branch.py --input examples/protein_input.tsv --out results/protein
"""
from __future__ import annotations
import argparse
import io
import re
import sys
from pathlib import Path
from typing import Any
from urllib.parse import quote
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from common import (Client, Variant, add_common_arguments, read_inputs, resolve,
                    save_run, vep, write_json, write_tsv)

UNIPROT = "https://rest.uniprot.org/uniprotkb"
AFDB = "https://alphafold.ebi.ac.uk/api/prediction"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def alphamissense_annotation(t: dict) -> dict:
    """REST/JSON plugins may be nested; older fixtures may have flat fields."""
    for key in ("alphamissense", "AlphaMissense"):
        value = t.get(key)
        if isinstance(value, dict):
            return value
    return {k: t[k] for k in ("am_pathogenicity", "am_class", "am_transcript_id", "am_uniprot_id", "am_protein_variant") if k in t}


def transcript_rank(t: dict) -> tuple:
    return (not bool(t.get("mane_select")), not bool(t.get("canonical")), t.get("transcript_id", ""))


def preferred_transcripts(transcripts: list[dict]) -> set[str]:
    """One preferred MISSENSE transcript per gene; all rows are still exported."""
    groups: dict[str, list[dict]] = {}
    for t in transcripts:
        groups.setdefault(t.get("gene_id") or t.get("transcript_id", ""), []).append(t)
    return {min(group, key=transcript_rank).get("transcript_id", "") for group in groups.values()}


def accessions(t: dict) -> list[str]:
    out = []
    for key in ("swissprot", "trembl"):
        values = t.get(key, [])
        if isinstance(values, (str, dict)):
            values = [values]
        for value in values:
            if isinstance(value, dict):
                value = value.get("id") or value.get("primary_id") or value.get("accession", "")
            value = str(value)
            # Isoform accessions are kept as candidates. Canonical is also tried, but
            # neither will be accepted without an exact full sequence match.
            for token in re.split(r"[,;\s]+", value):
                token = re.sub(r"\.\d+$", "", token)  # optional UniProt sequence-version suffix
                if re.fullmatch(r"[A-Z0-9]{6,10}(?:-\d+)?", token):
                    for item in (token, token.split("-")[0]):
                        if item not in out:
                            out.append(item)
    return out


def check_reference_position(sequence: str, pos: int, ref_aa: str) -> None:
    if not 1 <= pos <= len(sequence) or sequence[pos - 1] != ref_aa:
        raise ValueError("Protein reference residue does not match VEP; no structure mapping made")


def features_at_position(data: dict, pos: int) -> list[dict]:
    features = []
    for f in data.get("features", []):
        location = f.get("location", {})
        start, end = location.get("start", {}), location.get("end", {})
        if start.get("modifier", "EXACT") != "EXACT" or end.get("modifier", "EXACT") != "EXACT":
            continue
        lo, hi = start.get("value"), end.get("value")
        if isinstance(lo, int) and isinstance(hi, int) and lo <= pos <= hi:
            features.append({"type": f.get("type"), "start": lo, "end": hi,
                             "description": f.get("description", ""), "evidences": f.get("evidences", [])})
    return features


def comments_by_type(data: dict, category: str) -> list[dict]:
    return [c for c in data.get("comments", []) if c.get("commentType") == category]


def protein_name(data: dict) -> str:
    desc = data.get("proteinDescription", {})
    recommended = desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if recommended:
        return recommended
    names = desc.get("submissionNames", [])
    return names[0].get("fullName", {}).get("value", "") if names else ""


def pdb_residue_context(pdb_text: str, expected_seq: str, pos: int) -> dict:
    """Map by full chain sequence, then return the actual PDB residue identifier.
    Never assume that transcript residue number equals PDB resSeq number.
    """
    structure = PDBParser(QUIET=True).get_structure("AF", io.StringIO(pdb_text))
    model = next(structure.get_models())
    for chain in model:
        residues = [r for r in chain.get_residues() if is_aa(r, standard=True)]
        seq = "".join(seq1(r.get_resname()) for r in residues)
        if seq != expected_seq:
            continue
        if not 1 <= pos <= len(residues):
            raise ValueError("Variant residue outside matched structure chain")
        residue = residues[pos - 1]
        atom = residue["CA"] if "CA" in residue else next(residue.get_atoms())
        return {"chain": chain.id, "pdb_resseq": int(residue.id[1]),
                "pdb_insertion_code": residue.id[2].strip(),
                "reference_amino_acid": seq[pos - 1],
                "plddt": float(atom.bfactor),
                "plddt_interpretation": "local structure confidence, not variant pathogenicity"}
    raise ValueError("No complete PDB chain exactly matches the verified reference protein")


def alphafold_context(client: Client, accession: str, seq: str, pos: int, out: Path) -> dict:
    api_url = f"{AFDB}/{quote(accession, safe='')}"
    records = client.json(api_url)
    if not isinstance(records, list):
        raise ValueError("Unexpected AlphaFold API response schema")
    attempts = []
    for entry in records:
        af_seq = entry.get("uniprotSequence") or entry.get("sequence") or ""
        start = entry.get("uniprotStart", entry.get("sequenceStart", 1))
        end = entry.get("uniprotEnd", entry.get("sequenceEnd", len(af_seq)))
        if entry.get("isComplex") or af_seq != seq or int(start) != 1 or int(end) != len(seq):
            continue
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            continue
        model_id = entry.get("entryId") or entry.get("modelEntityId") or accession
        try:
            pdb_text = client.text(pdb_url, accept="text/plain")
            context = pdb_residue_context(pdb_text, seq, pos)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", str(model_id))
            local = out / "structures" / (safe_name + ".pdb")
            local.parent.mkdir(parents=True, exist_ok=True)
            local.write_text(pdb_text, encoding="utf-8")
            return {"status": "mapped_exact_full_sequence", "model_id": model_id,
                    "model_version": entry.get("latestVersion"), "model_date": entry.get("modelCreatedDate"),
                    "pdb_url": pdb_url, "cif_url": entry.get("cifUrl"),
                    "viewer_url": f"https://alphafold.ebi.ac.uk/entry/{accession}",
                    "pdb_local": str(local), "metadata_url": api_url, **context}
        except Exception as exc:
            attempts.append({"model_id": model_id, "error": str(exc)})
    return {"status": "no_verified_complete_model", "attempts": attempts,
            "note": "Missing/fragmentary/mismatched models are not evidence of a benign variant"}


def enrich(client: Client, t: dict, out: Path) -> dict:
    peptide = t.get("protein_id")
    pos, end = t.get("protein_start"), t.get("protein_end")
    aas = str(t.get("amino_acids", "")).split("/")
    if not peptide or not isinstance(pos, int) or pos != end or len(aas) != 2:
        return {"status": "unsupported_protein_change", "note": "Require a single amino-acid replacement"}
    ref_aa, alt_aa = aas
    if ref_aa not in STANDARD_AA or alt_aa not in STANDARD_AA:
        return {"status": "unsupported_amino_acid"}
    seq_data = client.ensembl(f"/sequence/id/{quote(peptide, safe='')}", {"type": "protein"})
    seq = str(seq_data.get("seq", "")).upper().rstrip("*")
    check_reference_position(seq, pos, ref_aa)
    candidates, attempts = accessions(t), []
    for accession in candidates:
        url = f"{UNIPROT}/{quote(accession, safe='')}.json"
        try:
            data = client.json(url)
            uniprot_seq = str(data.get("sequence", {}).get("value", "")).upper()
            if uniprot_seq != seq:
                attempts.append({"accession": accession, "reason": "full_sequence_mismatch"})
                continue
            result = {"status": "uniprot_exact_full_sequence", "uniprot_accession": accession,
                      "uniprot_url": url, "protein_name": protein_name(data),
                      "protein_length": len(seq), "protein_position": pos,
                      "amino_acid_change": f"{ref_aa}{pos}{alt_aa}",
                      "function_annotations": comments_by_type(data, "FUNCTION"),
                      "protein_level_disease_annotations": comments_by_type(data, "DISEASE"),
                      "disease_annotation_scope": "protein-level background; NOT evidence for this allele",
                      "overlapping_features": features_at_position(data, pos), "attempts": attempts}
            try:
                result["structure"] = alphafold_context(client, accession, seq, pos, out)
            except Exception as exc:
                result["structure"] = {"status": "unavailable", "error": str(exc)}
            return result
        except Exception as exc:
            attempts.append({"accession": accession, "reason": str(exc)})
    return {"status": "no_exact_uniprot_sequence_match", "candidates": candidates, "attempts": attempts,
            "note": "No approximate alignment or automatic isoform position transfer was attempted"}


def run(args: argparse.Namespace) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    client = Client(args.cache, args.offline, args.refresh)
    rows, variants_summary, errors = [], [], []
    seen = set()
    for index, input_row in enumerate(read_inputs(args.input), 1):
        try:
            variants = resolve(client, input_row)
        except Exception as exc:
            errors.append({"input_row": index, "input": input_row, "stage": "resolve", "error": str(exc)})
            continue
        for v in variants:
            if v.key in seen:
                continue
            seen.add(v.key)
            print(f"Annotating {v.key} ({v.rsid or 'coordinate input'})", file=sys.stderr)
            try:
                records = vep(client, v, not args.no_alphamissense)
                raw_path = out / "vep" / f"{v.chrom}_{v.pos}_{v.ref}_{v.alt}.json"
                write_json(raw_path, records)
                transcripts = [t for record in records for t in record.get("transcript_consequences", [])
                               if t.get("variant_allele") == v.alt and "missense_variant" in t.get("consequence_terms", [])]
                selected = preferred_transcripts(transcripts)
                variants_summary.append({**v.row(), "status": "annotated" if transcripts else "no_missense_consequence",
                                         "n_missense_transcripts": len(transcripts),
                                         "most_severe_consequences": sorted({r.get("most_severe_consequence", "") for r in records}),
                                         "vep_json": str(raw_path)})
                for t in sorted(transcripts, key=transcript_rank):
                    preferred = t.get("transcript_id", "") in selected
                    am = alphamissense_annotation(t)
                    row = {**v.row(), "gene": t.get("gene_symbol"), "gene_id": t.get("gene_id"),
                           "transcript": t.get("transcript_id"), "mane_select": t.get("mane_select"),
                           "canonical": t.get("canonical", 0), "preferred_missense_transcript": preferred,
                           "protein_id": t.get("protein_id"), "protein_position": t.get("protein_start"),
                           "amino_acids": t.get("amino_acids"), "hgvsp": t.get("hgvsp"),
                           "consequences": t.get("consequence_terms"), "domains": t.get("domains", []),
                           "vep_am_pathogenicity": am.get("am_pathogenicity"), "vep_am_class": am.get("am_class"),
                           "am_score_scope": "VEP-returned genomic-allele prediction; exact AlphaMissense transcript not guaranteed",
                           "am_note": "Prediction, not ClinVar classification, phenotype or disease probability",
                           "am_status": "returned_by_VEP" if am.get("am_pathogenicity") is not None else "not_returned",
                           "am_source": "Ensembl VEP AlphaMissense plugin", "am_raw_annotation": am}
                    if getattr(args, "vep_only", False):
                        context = {"status": "skipped_vep_only", "note": "VEP and AlphaMissense only; no UniProt or structure"}
                    elif preferred or args.all_transcripts:
                        try:
                            context = enrich(client, t, out)
                        except Exception as exc:
                            context = {"status": "annotation_error", "error": str(exc)}
                            errors.append({"variant_key": v.key, "transcript": t.get("transcript_id"),
                                           "stage": "protein_enrichment", "error": str(exc)})
                    else:
                        context = {"status": "not_selected_for_deep_annotation"}
                    row["protein_context"] = context
                    structure = context.get("structure", {})
                    row.update({"uniprot": context.get("uniprot_accession"), "protein_name": context.get("protein_name"),
                                "protein_mapping_status": context.get("status"),
                                "structure_status": structure.get("status", "not_available"),
                                "plddt": structure.get("plddt"), "pdb_url": structure.get("pdb_url"),
                                "pdb_chain": structure.get("chain"), "pdb_resseq": structure.get("pdb_resseq"),
                                "overlapping_features": context.get("overlapping_features", [])})
                    rows.append(row)
            except Exception as exc:
                variants_summary.append({**v.row(), "status": "annotation_error", "error": str(exc)})
                errors.append({"variant_key": v.key, "stage": "vep", "error": str(exc)})
    write_json(out / "protein_results.json", rows)
    write_json(out / "variant_status.json", variants_summary)
    write_tsv(out / "protein_results.tsv", rows, ["variant_key", "rsid", "chrom", "pos", "ref", "alt", "gene",
        "transcript", "mane_select", "canonical", "preferred_missense_transcript", "protein_id", "protein_position",
        "amino_acids", "hgvsp", "vep_am_pathogenicity", "vep_am_class", "am_status", "am_score_scope", "uniprot", "protein_name",
        "protein_mapping_status", "structure_status", "plddt", "pdb_chain", "pdb_resseq", "pdb_url", "overlapping_features"])
    save_run(out, client, args, errors, {"n_resolved_variants": len(seen), "n_missense_rows": len(rows),
        "n_errors": len(errors), "note": "Research prototype; no mutant structure or phenotype prediction"})
    print(f"Saved {len(rows)} missense transcript rows; {len(errors)} logged errors to {out}")
    return 0 if any(s.get("status") in {"annotated", "no_missense_consequence"} for s in variants_summary) else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument("--all-transcripts", action="store_true", help="Enrich every missense transcript, not just one per gene")
    parser.add_argument("--vep-only", action="store_true", help="Annotate with VEP without downloading proteins or structures")
    parser.add_argument("--no-alphamissense", action="store_true", help="Disable optional VEP AlphaMissense plugin explicitly")
    args = parser.parse_args()
    try:
        return run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
