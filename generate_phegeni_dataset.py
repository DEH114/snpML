"""Build an NCBI PheGenI association dataset for an rsID panel.

The official PheGenI bulk association table combines NHGRI GWAS Catalog
records with selected top results from dbGaP analyses. Matching is performed
only against the complete numeric value in PheGenI's ``SNP rs`` column and
against both the supplied and current dbSNP identifiers. Every matched source
row is retained losslessly alongside an analysis-friendly normalized table.
By default the exact QIAGEN 20-rsID panel is required; use
--allow-custom-panel for another validated, unique rsID list.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any, Iterable

import requests

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import Client, utc_now, write_json, write_tsv
from generate_gwas_catalog_dataset import (
    RESOLUTION_FIELDS,
    resolve_dbsnp_ids,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "phegeni_dataset"
PHEGENI_TABLE_URL = (
    "https://www.ncbi.nlm.nih.gov/projects/gap/eqtl/"
    "EpiViewBE.cgi?type=dl.tab"
)
TABLE_FILENAME = "PheGenI_Association_full.tab"
USER_AGENT = "BioHack-SNP-Toolkit/1.0 (research prototype)"
DEFAULT_MAX_DOWNLOAD_BYTES = 100 * 1024 * 1024

SOURCE_FIELDS = [
    "#", "Trait", "SNP rs", "Context", "Gene", "Gene ID", "Gene 2",
    "Gene ID 2", "Chromosome", "Location", "P-Value", "Source",
    "PubMed", "Analysis ID", "Study ID", "Study Name",
]
REQUIRED_SOURCE_FIELDS = set(SOURCE_FIELDS)

ASSOCIATION_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "source_row_number",
    "source_record_sha256", "match_basis", "matched_rsid",
    "phegeni_row_id", "trait", "source_snp_rsid", "source_snp_numeric",
    "context", "gene", "gene_id", "gene_2", "gene_id_2", "chromosome",
    "location", "p_value", "evidence_source", "pubmed_id", "analysis_id",
    "study_id", "study_name", "source_url",
]

STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "dbsnp_status", "status",
    "n_associations", "n_nhgri_associations", "n_dbgap_associations",
    "n_other_source_associations", "n_unique_traits", "n_unique_genes",
    "n_unique_publications", "n_unique_analyses", "n_unique_studies",
    "matched_source_rsids", "error",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_record_digest(row: dict[str, str]) -> str:
    encoded = json.dumps(
        row, ensure_ascii=False, separators=(",", ":"), sort_keys=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def exact_numeric_rsid(value: Any) -> str | None:
    """Return an rsID only when *value* is one complete numeric identifier."""
    text = str(value or "").strip()
    if not re.fullmatch(r"[1-9][0-9]*", text):
        return None
    return f"rs{text}"


def inspect_table(path: Path, max_bytes: int) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"PheGenI association table not found: {path}")
    size = path.stat().st_size
    if size < 1:
        raise ValueError("PheGenI association table is empty")
    if size > max_bytes:
        raise ValueError(
            f"PheGenI table is {size:,} bytes, exceeding the "
            f"{max_bytes:,}-byte safety limit"
        )
    with path.open("rb") as handle:
        prefix = handle.read(4096).lstrip().lower()
    if prefix.startswith((b"<!doctype html", b"<html", b"<?xml")):
        raise ValueError("PheGenI server returned HTML/XML instead of tabular data")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
    except UnicodeDecodeError as exc:
        raise ValueError("PheGenI table is not valid UTF-8 text") from exc
    if len(set(header)) != len(header):
        raise ValueError("PheGenI association table has duplicate headers")
    missing = sorted(REQUIRED_SOURCE_FIELDS - set(header))
    if missing:
        raise ValueError(f"PheGenI association table is missing columns: {missing}")
    return {
        "table_path": str(path.resolve()),
        "table_bytes": size,
        "table_sha256": sha256_file(path),
        "source_header": header,
    }


def download_table(url: str, destination: Path, max_bytes: int) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    digest = hashlib.sha256()
    total = 0
    try:
        with requests.get(
            url,
            stream=True,
            timeout=(15, 180),
            headers={"User-Agent": USER_AGENT},
        ) as response:
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "")
            if "text/html" in content_type.lower():
                raise ValueError("PheGenI server returned HTML instead of tabular data")
            declared = response.headers.get("Content-Length", "")
            if declared.isdigit() and int(declared) > max_bytes:
                raise ValueError(
                    f"PheGenI download declares {int(declared):,} bytes, exceeding "
                    f"the {max_bytes:,}-byte safety limit"
                )
            with partial.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > max_bytes:
                        raise ValueError(
                            "PheGenI download exceeded the configured size limit"
                        )
                    digest.update(chunk)
                    handle.write(chunk)
            provenance = {
                "requested_url": url,
                "resolved_url": response.url,
                "retrieved_at_utc": utc_now(),
                "downloaded_bytes": total,
                "download_sha256": digest.hexdigest(),
                "content_type": content_type,
                "content_disposition": response.headers.get("Content-Disposition"),
                "last_modified": response.headers.get("Last-Modified"),
                "etag": response.headers.get("ETag"),
            }
        inspect_table(partial, max_bytes)
        partial.replace(destination)
        return provenance
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def acquire_table(args: argparse.Namespace) -> tuple[Path, dict]:
    local = getattr(args, "table", None)
    if local:
        path = Path(local)
        provenance: dict[str, Any] = {
            "source_mode": "local",
            "requested_url": "",
        }
    else:
        path = Path(args.cache) / "phegeni" / TABLE_FILENAME
        sidecar = path.with_suffix(path.suffix + ".source.json")
        if args.refresh or not path.exists():
            if args.offline:
                raise FileNotFoundError(
                    f"Offline PheGenI cache miss: {path}. Supply --table or run "
                    "once without --offline."
                )
            provenance = download_table(
                PHEGENI_TABLE_URL, path, args.max_download_bytes,
            )
            provenance["source_mode"] = "downloaded"
            write_json(sidecar, provenance)
        else:
            provenance = {
                "source_mode": "cache",
                "requested_url": PHEGENI_TABLE_URL,
            }
            if sidecar.is_file():
                cached = json.loads(sidecar.read_text(encoding="utf-8"))
                if isinstance(cached, dict):
                    provenance.update(cached)
                    provenance["source_mode"] = "cache"
    inspected = inspect_table(path, args.max_download_bytes)
    provenance.update(inspected)
    return path, provenance


def normalized_association(
    panel_rsid: str,
    resolution: dict,
    rank: int,
    source_row_number: int,
    digest: str,
    row: dict[str, str],
    source_rsid: str,
) -> dict:
    current = resolution["current_rsid"]
    if source_rsid == current and current != panel_rsid:
        match_basis = "current_rsid"
    else:
        match_basis = "input_rsid"
    return {
        "rank": rank,
        "rsid": panel_rsid,
        "current_rsid": current,
        "was_merged": current != panel_rsid,
        "source_row_number": source_row_number,
        "source_record_sha256": digest,
        "match_basis": match_basis,
        "matched_rsid": source_rsid,
        "phegeni_row_id": row["#"],
        "trait": row["Trait"],
        "source_snp_rsid": source_rsid,
        "source_snp_numeric": row["SNP rs"].strip(),
        "context": row["Context"],
        "gene": row["Gene"],
        "gene_id": row["Gene ID"],
        "gene_2": row["Gene 2"],
        "gene_id_2": row["Gene ID 2"],
        "chromosome": row["Chromosome"],
        "location": row["Location"],
        "p_value": row["P-Value"],
        "evidence_source": row["Source"],
        "pubmed_id": row["PubMed"],
        "analysis_id": row["Analysis ID"],
        "study_id": row["Study ID"],
        "study_name": row["Study Name"],
        "source_url": (
            "https://www.ncbi.nlm.nih.gov/gap/phegeni?rs="
            f"{source_rsid.removeprefix('rs')}"
        ),
    }


def parse_table(
    path: Path,
    rsids: list[str],
    resolutions: list[dict],
) -> tuple[list[dict], list[dict], list[str], int]:
    resolution_by_input = {row["input_rsid"]: row for row in resolutions}
    aliases: dict[str, list[str]] = {}
    rank_by_rsid = {rsid: rank for rank, rsid in enumerate(rsids, 1)}
    for rsid in rsids:
        current = resolution_by_input[rsid]["current_rsid"]
        for alias in dict.fromkeys([rsid, current]):
            aliases.setdefault(alias, []).append(rsid)

    associations: list[dict] = []
    matched_source_rows: list[dict] = []
    scanned = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if len(set(header)) != len(header):
            raise ValueError("PheGenI association table has duplicate headers")
        missing = sorted(REQUIRED_SOURCE_FIELDS - set(header))
        if missing:
            raise ValueError(f"PheGenI association table is missing columns: {missing}")
        for source_row_number, raw in enumerate(reader, 2):
            scanned += 1
            if None in raw or any(raw.get(name) is None for name in header):
                raise ValueError(
                    f"Malformed PheGenI row {source_row_number}: wrong field count"
                )
            row = {name: raw.get(name, "") or "" for name in header}
            source_rsid = exact_numeric_rsid(row["SNP rs"])
            if source_rsid is None:
                raise ValueError(
                    f"Malformed PheGenI row {source_row_number}: SNP rs must be "
                    f"one exact positive integer, got {row['SNP rs']!r}"
                )
            matching_inputs = aliases.get(source_rsid, [])
            if not matching_inputs:
                continue
            ordered_inputs = sorted(matching_inputs, key=rank_by_rsid.__getitem__)
            digest = source_record_digest(row)
            matched_source_rows.append({
                "source_row_number": source_row_number,
                "source_record_sha256": digest,
                "matched_panel_rsids": ordered_inputs,
                "record": row,
            })
            for panel_rsid in ordered_inputs:
                associations.append(normalized_association(
                    panel_rsid,
                    resolution_by_input[panel_rsid],
                    rank_by_rsid[panel_rsid],
                    source_row_number,
                    digest,
                    row,
                    source_rsid,
                ))
    associations.sort(key=lambda row: (row["rank"], row["source_row_number"]))
    return associations, matched_source_rows, header, scanned


def build_statuses(
    rsids: list[str], resolutions: list[dict], associations: list[dict],
) -> list[dict]:
    resolution_by_input = {row["input_rsid"]: row for row in resolutions}
    by_rsid = {rsid: [] for rsid in rsids}
    for association in associations:
        by_rsid[association["rsid"]].append(association)
    statuses = []
    for rank, rsid in enumerate(rsids, 1):
        rows = by_rsid[rsid]
        resolution = resolution_by_input[rsid]
        nhgri = sum(row["evidence_source"].casefold() == "nhgri" for row in rows)
        dbgap = sum(row["evidence_source"].casefold() == "dbgap" for row in rows)
        statuses.append({
            "rank": rank,
            "rsid": rsid,
            "current_rsid": resolution["current_rsid"],
            "was_merged": resolution["current_rsid"] != rsid,
            "dbsnp_status": resolution.get("status", ""),
            "status": "associations_found" if rows else "no_phegeni_association",
            "n_associations": len(rows),
            "n_nhgri_associations": nhgri,
            "n_dbgap_associations": dbgap,
            "n_other_source_associations": len(rows) - nhgri - dbgap,
            "n_unique_traits": len({row["trait"] for row in rows if row["trait"]}),
            "n_unique_genes": len({
                gene
                for row in rows
                for gene in (row["gene"], row["gene_2"])
                if gene
            }),
            "n_unique_publications": len({
                row["pubmed_id"] for row in rows if row["pubmed_id"] not in {"", "0"}
            }),
            "n_unique_analyses": len({
                row["analysis_id"] for row in rows if row["analysis_id"] not in {"", "0"}
            }),
            "n_unique_studies": len({
                row["study_id"] for row in rows if row["study_id"] not in {"", "0"}
            }),
            "matched_source_rsids": sorted({row["source_snp_rsid"] for row in rows}),
            "error": "",
        })
    return statuses


def sqlite_value(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, bool):
        return int(value)
    return value


def create_table(
    database: sqlite3.Connection,
    name: str,
    fields: list[str],
    integer_fields: set[str],
    primary_key: str = "",
) -> None:
    columns = ", ".join(
        f'"{field}" {"INTEGER" if field in integer_fields else "TEXT"}'
        for field in fields
    )
    suffix = f", PRIMARY KEY ({primary_key})" if primary_key else ""
    database.execute(f'CREATE TABLE "{name}" ({columns}{suffix})')


def insert_rows(
    database: sqlite3.Connection,
    name: str,
    fields: list[str],
    rows: Iterable[dict],
) -> None:
    columns = ", ".join(f'"{field}"' for field in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{name}" ({columns}) VALUES ({placeholders})',
        [tuple(sqlite_value(row.get(field, "")) for field in fields) for row in rows],
    )


def write_sqlite(
    path: Path,
    resolutions: list[dict],
    associations: list[dict],
    raw_rows: list[dict],
    statuses: list[dict],
    manifest: dict,
) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    raw_fields = [
        "source_row_number", "source_record_sha256", "matched_panel_rsids", "record",
    ]
    try:
        with database:
            create_table(
                database, "dbsnp_resolution", RESOLUTION_FIELDS,
                {"rank", "was_merged"}, '"input_rsid"',
            )
            create_table(
                database, "phegeni_associations", ASSOCIATION_FIELDS,
                {"rank", "was_merged", "source_row_number"},
                '"rsid", "source_row_number"',
            )
            create_table(
                database, "matched_source_records", raw_fields,
                {"source_row_number"}, '"source_row_number"',
            )
            status_integers = {"rank", "was_merged"} | {
                field for field in STATUS_FIELDS if field.startswith("n_")
            }
            create_table(
                database, "panel_status", STATUS_FIELDS, status_integers, '"rsid"',
            )
            database.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)"
            )
            insert_rows(database, "dbsnp_resolution", RESOLUTION_FIELDS, resolutions)
            insert_rows(
                database, "phegeni_associations", ASSOCIATION_FIELDS, associations,
            )
            insert_rows(database, "matched_source_records", raw_fields, raw_rows)
            insert_rows(database, "panel_status", STATUS_FIELDS, statuses)
            database.execute(
                "INSERT INTO metadata VALUES (?, ?)",
                ("manifest", json.dumps(manifest, ensure_ascii=False, sort_keys=True)),
            )
            integrity = database.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise ValueError(f"SQLite integrity check failed: {integrity}")
    except Exception:
        database.close()
        partial.unlink(missing_ok=True)
        raise
    else:
        database.close()
    partial.replace(path)


def write_readme(out: Path, manifest: dict) -> None:
    text = f"""# PheGenI dataset for a {manifest['n_input_rsids']}-rsID panel

This folder is a reproducible, public association-data snapshot from NCBI's
Phenotype-Genotype Integrator (PheGenI). It contains no patient-level data.

- Source table SHA-256: {manifest['source_provenance']['table_sha256']}
- Input rsIDs: {manifest['n_input_rsids']}
- rsIDs with one or more associations: {manifest['n_rsids_with_associations']}
- Normalized rsID-association relationships: {manifest['n_association_relations']}
- Losslessly retained source rows: {manifest['n_matched_source_rows']}
- NHGRI-source relationships: {manifest['n_nhgri_association_relations']}
- dbGaP-source relationships: {manifest['n_dbgap_association_relations']}

## Files

`phegeni_associations.tsv` contains one normalized relationship between a panel
rsID and a matching PheGenI row. `matched_source_records.json` retains every
matching source row, all original columns, and all original strings without
numeric coercion. `panel_status.tsv` includes every input in original order,
including rsIDs with no match. `dbsnp_resolution.tsv` records current IDs used
for matching. `source_provenance.json` and `run_manifest.json` record hashes,
retrieval details, and counts. `phegeni.sqlite` mirrors the normalized tables.

Matching is exact: the complete numeric `SNP rs` value must equal the numeric
part of the supplied rsID or its current dbSNP ID. Substring, list, range, and
free-text matches are never inferred.

PheGenI is an integrator, not an independent seventh experiment. Its table
combines NHGRI GWAS Catalog records with selected top results from dbGaP, so
some NHGRI records can overlap the separate GWAS Catalog dataset. Do not add
counts across those datasets as if they were independent evidence. A reported
association is not proof of causality, a diagnosis, or an individual risk
estimate. PheGenI's bulk file does not expose an explicit release identifier;
the recorded SHA-256 is the authoritative snapshot identifier.

Official source:
https://www.ncbi.nlm.nih.gov/gap/phegeni
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def run(args: argparse.Namespace, client: Client | None = None) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    if args.max_download_bytes < 1:
        raise ValueError("--max-download-bytes must be positive")
    rsids = select_panel_rsids(
        args.input,
        require_exact=not getattr(args, "allow_custom_panel", False),
    )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for stale_name in ("phegeni.sqlite", "phegeni.sqlite.part"):
        (out / stale_name).unlink(missing_ok=True)
    (out / "raw" / "dbsnp_esummary.json").unlink(missing_ok=True)

    warnings: list[dict] = []
    client = client or Client(args.cache, args.offline, args.refresh)
    source_path, source_provenance = acquire_table(args)
    resolutions, resolution_source, dbsnp_raw = resolve_dbsnp_ids(
        args, rsids, client, warnings,
    )
    if dbsnp_raw is not None:
        write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)

    associations, raw_rows, source_header, scanned = parse_table(
        source_path, rsids, resolutions,
    )
    statuses = build_statuses(rsids, resolutions, associations)
    input_path = Path(args.input)
    nhgri_count = sum(
        row["evidence_source"].casefold() == "nhgri" for row in associations
    )
    dbgap_count = sum(
        row["evidence_source"].casefold() == "dbgap" for row in associations
    )
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source_database": "NCBI Phenotype-Genotype Integrator (PheGenI)",
        "source_table_url": source_provenance.get("requested_url", ""),
        "source_release": "not published by the bulk endpoint; use table SHA-256",
        "dbsnp_resolution_source": resolution_source,
        "n_input_rsids": len(rsids),
        "n_merged_input_rsids": sum(bool(row["was_merged"]) for row in resolutions),
        "n_source_rows_scanned": scanned,
        "n_matched_source_rows": len(raw_rows),
        "n_association_relations": len(associations),
        "n_nhgri_association_relations": nhgri_count,
        "n_dbgap_association_relations": dbgap_count,
        "n_other_source_association_relations": (
            len(associations) - nhgri_count - dbgap_count
        ),
        "n_rsids_with_associations": sum(
            row["n_associations"] > 0 for row in statuses
        ),
        "n_warnings": len(warnings),
        "source_header": source_header,
        "source_provenance": source_provenance,
        "http_provenance": list(client.manifest.values()),
        "notes": [
            "SNP rs is matched only as one complete positive integer",
            "Both supplied and current dbSNP identifiers are considered",
            "All original source strings remain lossless in matched_source_records",
            "PheGenI integrates NHGRI GWAS Catalog and selected dbGaP top results",
            "PheGenI may overlap the standalone GWAS Catalog dataset",
            "Association does not establish causality or individual clinical risk",
        ],
    }

    raw_tsv_rows = [{
        "source_row_number": item["source_row_number"],
        "source_record_sha256": item["source_record_sha256"],
        "matched_panel_rsids": item["matched_panel_rsids"],
        **item["record"],
    } for item in raw_rows]
    raw_tsv_fields = [
        "source_row_number", "source_record_sha256", "matched_panel_rsids",
        *source_header,
    ]

    write_json(out / "source_provenance.json", source_provenance)
    write_json(out / "dbsnp_resolution.json", resolutions)
    write_tsv(out / "dbsnp_resolution.tsv", resolutions, RESOLUTION_FIELDS)
    write_json(out / "matched_source_records.json", raw_rows)
    write_tsv(out / "matched_source_records.tsv", raw_tsv_rows, raw_tsv_fields)
    write_json(out / "phegeni_associations.json", associations)
    write_tsv(out / "phegeni_associations.tsv", associations, ASSOCIATION_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "errors.json", warnings)
    write_json(out / "run_manifest.json", manifest)
    write_sqlite(
        out / "phegeni.sqlite", resolutions, associations, raw_rows, statuses, manifest,
    )
    write_readme(out, manifest)
    print(
        f"PheGenI dataset complete: {len(rsids)} input rsIDs, "
        f"{len(associations)} association relationships from "
        f"{len(raw_rows)} source rows. See {out}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", default=str(DEFAULT_INPUT),
        help="rsID-only TSV or one-rsID-per-line input",
    )
    parser.add_argument(
        "--allow-custom-panel", action="store_true",
        help="Accept any validated, unique rsID list instead of the exact QIAGEN panel",
    )
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument(
        "--table", help="Local official PheGenI Association_full.tab (no download)",
    )
    parser.add_argument("--cache", default=".cache", help="Download/API cache directory")
    parser.add_argument(
        "--offline", action="store_true",
        help="Use --table or the cached PheGenI table and cached API data only",
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Refresh the official table and API cache",
    )
    parser.add_argument(
        "--dbsnp-resolution",
        help=("Existing dbsnp_resolution.json; by default reuse the sibling "
              "ClinVar dataset, then try NCBI ESummary"),
    )
    parser.add_argument(
        "--max-download-bytes", type=int, default=DEFAULT_MAX_DOWNLOAD_BYTES,
        help="Maximum accepted downloaded, cached, or local table size",
    )
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
