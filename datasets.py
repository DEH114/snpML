"""Download public resources and construct small, reproducible SNP demo panels.
Only downloads public, variant-level data, never patient VCFs.
ClinVar classes are source annotations for retrieval/demo, not independent ground truth.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
import re
import sys
import zipfile
import requests
from common import CHROMS, DNA, utc_now, write_json, write_tsv

SOURCES = {
    "clinvar": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261101.vcf.gz",
    "gwas": "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip",
}
INPUT_FIELDS = ["rsid", "chrom", "pos", "ref", "alt", "assembly"]


def download(source: str, destination: Path) -> None:
    url = SOURCES[source]
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    digest = hashlib.sha256()
    size = 0
    try:
        with requests.get(url, stream=True, timeout=(15, 120),
                          headers={"User-Agent": "BioHack-SNP-Toolkit/1.0"}) as response:
            response.raise_for_status()
            if "text/html" in response.headers.get("Content-Type", "").lower():
                raise ValueError("Server returned HTML instead of the requested dataset")
            with partial.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if not chunk:
                        continue
                    digest.update(chunk)
                    size += len(chunk)
                    handle.write(chunk)
            meta = {"source": source, "requested_url": url, "resolved_url": response.url,
                    "retrieved_at_utc": utc_now(), "bytes": size, "sha256": digest.hexdigest(),
                    "last_modified": response.headers.get("Last-Modified")}
        # Basic archive integrity check before declaring a successful download.
        if destination.suffix == ".zip":
            with zipfile.ZipFile(partial) as archive:
                if archive.testzip() is not None:
                    raise ValueError("ZIP failed CRC check")
        elif destination.suffix == ".gz":
            with gzip.open(partial, "rb") as handle:
                handle.read(1024)
        partial.replace(destination)
        write_json(str(destination) + ".source.json", meta)
        print(f"Downloaded {size:,} bytes to {destination}")
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def parse_info(text: str) -> dict[str, str]:
    result = {}
    for token in text.split(";"):
        key, sep, value = token.partition("=")
        result[key] = value if sep else ""
    return result


def clinvar_class(info: dict) -> str | None:
    """Conservative >=2-star germline classes, with no conflicting-interpretation records."""
    accepted_reviews = {"criteria_provided,_multiple_submitters,_no_conflicts",
                        "reviewed_by_expert_panel", "practice_guideline"}
    if info.get("CLNREVSTAT") not in accepted_reviews or info.get("CLNSIGCONF"):
        return None
    labels = set(re.split(r"[/|]", info.get("CLNSIG", "")))
    if labels and labels <= {"Pathogenic", "Likely_pathogenic"}:
        return "pathogenic_or_likely_pathogenic"
    if labels and labels <= {"Benign", "Likely_benign"}:
        return "benign_or_likely_benign"
    return None


def clinvar_panel(vcf: Path, out: Path, n: int, seed: int, missense_only: bool) -> None:
    if n < 1:
        raise ValueError("n-per-class must be positive")
    rng = random.Random(seed)
    labels = ["pathogenic_or_likely_pathogenic", "benign_or_likely_benign"]
    pools = {label: [] for label in labels}
    counts = {label: 0 for label in labels}
    opener = gzip.open if vcf.suffix == ".gz" else open
    confirmed_build = False
    saw_variant = False
    with opener(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                if line.lower().startswith(("##reference=", "##assembly=")) and "GRCh38" in line:
                    confirmed_build = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            if not confirmed_build:
                raise ValueError("VCF lacks a GRCh38 reference/assembly header; do not assume the genome build")
            saw_variant = True
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Malformed VCF record")
            chrom, pos, cv_id, ref, alt = fields[:5]
            chrom = chrom.removeprefix("chr")
            # Skip multiallelic records instead of incorrectly assigning aggregate CLNSIG to an ALT.
            if chrom not in CHROMS or ref not in DNA or alt not in DNA or ref == alt:
                continue
            info = parse_info(fields[7])
            label = clinvar_class(info)
            if label is None:
                continue
            if missense_only and "missense_variant" not in info.get("MC", "") and "SO:0001583" not in info.get("MC", ""):
                continue
            rs_numbers = [x for x in info.get("RS", "").split("|") if x.isdigit()]
            row = {"rsid": "rs" + rs_numbers[0] if len(rs_numbers) == 1 else "", "chrom": chrom,
                   "pos": int(pos), "ref": ref, "alt": alt, "assembly": "GRCh38", "source_class": label,
                   "clinvar_variation_id": cv_id, "clinical_significance": info.get("CLNSIG", ""),
                   "review_status": info.get("CLNREVSTAT", ""), "conditions": info.get("CLNDN", ""),
                   "molecular_consequence": info.get("MC", ""), "gene_info": info.get("GENEINFO", ""),
                   "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/"}
            # Reservoir sample across the ENTIRE file, not the first n rows/first chromosome.
            counts[label] += 1
            pool = pools[label]
            if len(pool) < n:
                pool.append(row)
            else:
                slot = rng.randrange(counts[label])
                if slot < n:
                    pool[slot] = row
    if not saw_variant:
        raise ValueError("VCF contained no variants")
    selected = [r for label in labels for r in pools[label]]
    if not selected:
        raise ValueError("No variants met the requested filters")
    write_tsv(out, selected, INPUT_FIELDS + ["source_class", "clinvar_variation_id", "clinical_significance",
        "review_status", "conditions", "molecular_consequence", "gene_info", "source_url"])
    write_json(str(out) + ".selection.json", {"input_vcf": str(vcf), "seed": seed, "n_requested_per_class": n,
        "eligible_counts": counts, "selected_counts": {k: len(v) for k, v in pools.items()},
        "missense_only": missense_only, "assembly": "GRCh38", "method": "full-file per-class reservoir sampling",
        "warning": "Demo/retrieval panel, not an independent predictor benchmark; preserve original variant-condition evidence"})
    print(f"Saved {len(selected)} variants to {out}; eligible counts: {counts}")


@contextmanager
def open_tsv(path: Path):
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            entries = [n for n in archive.namelist() if not n.endswith("/") and n.lower().endswith((".tsv", ".txt"))]
            if len(entries) != 1:
                raise ValueError("ZIP must contain one TSV/TXT; unpack and specify the intended table otherwise")
            with archive.open(entries[0]) as binary:
                with io.TextIOWrapper(binary, encoding="utf-8-sig") as handle:
                    yield handle
    elif path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
            yield handle
    else:
        with path.open(encoding="utf-8-sig") as handle:
            yield handle


def gwas_panel(source: Path, out: Path, trait_regex: str, n: int, seed: int) -> None:
    if n < 1:
        raise ValueError("n must be positive")
    pattern = re.compile(trait_regex, re.I)
    groups: dict[str, list[dict]] = {}
    with open_tsv(source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"SNPS", "DISEASE/TRAIT"} <= set(reader.fieldnames or []):
            raise ValueError("Expected GWAS Catalog all-associations TSV with SNPS and DISEASE/TRAIT headers")
        for row in reader:
            if not pattern.search(row.get("DISEASE/TRAIT", "")):
                continue
            snps = row.get("SNPS", "").strip()
            # Do not explode haplotype or SNP-by-SNP interaction associations into single-SNP claims.
            if not re.fullmatch(r"rs\d+", snps):
                continue
            groups.setdefault(snps, []).append(row)
    if not groups:
        raise ValueError("No single-rsID associations matched that trait regex")
    ids = random.Random(seed).sample(sorted(groups), min(n, len(groups)))
    write_tsv(out, [{"rsid": rsid} for rsid in ids], ["rsid"])
    evidence = {rsid: groups[rsid] for rsid in ids}
    write_json(str(out) + ".gwas_records.json", evidence)
    write_json(str(out) + ".selection.json", {"input": str(source), "trait_regex": trait_regex, "seed": seed,
        "unique_rsids_eligible": len(groups), "selected": len(ids),
        "warning": "rsID-level input only; risk alleles are retained in raw records, not assumed to be genomic ALT. Resolve and harmonize before making allele-specific trait claims. A GWAS hit need not be the causal variant."})
    print(f"Saved {len(ids)} rsIDs to {out}; original associations saved separately")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    dl = sub.add_parser("download")
    dl.add_argument("--source", choices=list(SOURCES), required=True)
    dl.add_argument("--out", required=True)
    cv = sub.add_parser("clinvar")
    cv.add_argument("--vcf", required=True)
    cv.add_argument("--out", required=True)
    cv.add_argument("--n-per-class", type=int, default=25)
    cv.add_argument("--seed", type=int, default=42)
    cv.add_argument("--missense-only", action="store_true")
    gw = sub.add_parser("gwas")
    gw.add_argument("--input", required=True)
    gw.add_argument("--out", required=True)
    gw.add_argument("--trait", default="LDL|low density lipoprotein", help="Case-insensitive regular expression")
    gw.add_argument("--n", type=int, default=20)
    gw.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    try:
        if args.command == "download":
            download(args.source, Path(args.out))
        elif args.command == "clinvar":
            clinvar_panel(Path(args.vcf), Path(args.out), args.n_per_class, args.seed, args.missense_only)
        else:
            gwas_panel(Path(args.input), Path(args.out), args.trait, args.n, args.seed)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
