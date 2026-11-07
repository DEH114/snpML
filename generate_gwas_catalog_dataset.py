"""Build a GWAS Catalog dataset for the exact 20 QIAGEN challenge rsIDs.

The official ontology-annotated association ZIP is streamed directly from the
archive (never extracted).  Matching is by exact rsID token against both the
supplied identifier and its current dbSNP identifier.  Every matched source
row is retained losslessly, while a second table records one panel-rsID to
association relationship per match.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sqlite3
import sys
import zipfile
from pathlib import Path
from typing import Any, Iterable

import requests

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import Client, utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "gwas_catalog_dataset"
DEFAULT_DBSNP_RESOLUTION = (
    DEFAULT_OUT.parent / "clinvar_dataset" / "dbsnp_resolution.json"
)
GWAS_ARCHIVE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip"
)
GWAS_METADATA_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/metadata"
ARCHIVE_FILENAME = "gwas-catalog-associations_ontology-annotated-full.zip"
USER_AGENT = "BioHack-SNP-Toolkit/1.0 (research prototype)"
RSID_TOKEN = re.compile(r"(?<![A-Za-z0-9])rs([0-9]+)(?![A-Za-z0-9])", re.I)

RESOLUTION_FIELDS = [
    "rank", "input_rsid", "current_rsid", "was_merged", "dbsnp_uid",
    "snp_class", "spdi", "chromosome", "chrpos_zero_based", "merge_text",
    "status", "error",
]
ASSOCIATION_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "source_row_number",
    "source_record_sha256", "match_basis", "matched_rsid",
    "matched_identifiers", "reported_rsids", "current_rsids",
    "risk_allele_rsids", "association_scope", "is_multi_variant",
    "strongest_snp_risk_allele", "date_added_to_catalog", "pubmed_id",
    "first_author", "publication_date", "journal", "publication_link",
    "study_title", "reported_trait", "initial_sample_description",
    "replication_sample_description", "region", "chromosome", "position",
    "reported_genes", "mapped_genes", "context", "intergenic",
    "risk_allele_frequency", "p_value", "pvalue_mlog", "pvalue_text",
    "or_or_beta", "confidence_interval", "platform", "cnv",
    "mapped_trait", "mapped_trait_uri", "study_accession",
    "genotyping_technology", "source_url",
]
STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "dbsnp_status", "status",
    "n_associations", "n_single_variant_associations",
    "n_multi_variant_associations", "n_matches_by_reported_rsid",
    "n_matches_by_current_rsid", "n_unique_studies", "n_unique_traits",
    "error",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_rsid_tokens(value: Any) -> list[str]:
    """Return distinct, normalized rsIDs without substring matches."""
    result: list[str] = []
    for match in RSID_TOKEN.finditer(str(value or "")):
        rsid = f"rs{match.group(1)}"
        if rsid not in result:
            result.append(rsid)
    return result


def current_rsid_tokens(value: Any) -> list[str]:
    """Parse a SNP_ID_CURRENT field, whose official files may omit ``rs``."""
    text = str(value or "").strip()
    result = exact_rsid_tokens(text)
    for numeric in re.findall(r"(?<![0-9])([0-9]+)(?![0-9])", text):
        rsid = f"rs{numeric}"
        if rsid not in result:
            result.append(rsid)
    return result


def source_record_digest(row: dict[str, str]) -> str:
    encoded = json.dumps(
        row, ensure_ascii=False, separators=(",", ":"), sort_keys=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def choose_archive_member(archive: zipfile.ZipFile,
                          max_uncompressed_bytes: int) -> zipfile.ZipInfo:
    total_uncompressed = sum(info.file_size for info in archive.infolist()
                             if not info.is_dir())
    if total_uncompressed > max_uncompressed_bytes:
        raise ValueError(
            f"GWAS Catalog ZIP expands to {total_uncompressed:,} bytes, exceeding "
            f"the {max_uncompressed_bytes:,}-byte safety limit"
        )
    candidates = [
        info for info in archive.infolist()
        if not info.is_dir()
        and not info.filename.replace("\\", "/").startswith("__MACOSX/")
        and Path(info.filename).suffix.lower() in {".tsv", ".txt"}
    ]
    if len(candidates) != 1:
        raise ValueError(
            "GWAS Catalog ZIP must contain exactly one TSV/TXT data member; "
            f"found {[item.filename for item in candidates]}"
        )
    member = candidates[0]
    if member.file_size < 1:
        raise ValueError("GWAS Catalog archive member is empty")
    if member.file_size > max_uncompressed_bytes:
        raise ValueError(
            f"GWAS Catalog member is {member.file_size:,} bytes, exceeding the "
            f"{max_uncompressed_bytes:,}-byte safety limit"
        )
    return member


def inspect_archive(path: Path, max_uncompressed_bytes: int) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"GWAS Catalog archive not found: {path}")
    try:
        with zipfile.ZipFile(path) as archive:
            member = choose_archive_member(archive, max_uncompressed_bytes)
            failed_member = archive.testzip()
            if failed_member is not None:
                raise ValueError(f"GWAS Catalog ZIP failed CRC check: {failed_member}")
            member_info = {
                "member_name": member.filename,
                "member_uncompressed_bytes": member.file_size,
                "member_compressed_bytes": member.compress_size,
                "member_crc32": f"{member.CRC:08x}",
            }
    except zipfile.BadZipFile as exc:
        raise ValueError(f"Invalid GWAS Catalog ZIP: {path}") from exc
    return {
        "archive_path": str(path.resolve()),
        "archive_bytes": path.stat().st_size,
        "archive_sha256": sha256_file(path),
        "archive_crc_check": "ok",
        **member_info,
    }


def download_archive(url: str, destination: Path, max_download_bytes: int,
                     max_uncompressed_bytes: int) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    digest = hashlib.sha256()
    total = 0
    try:
        with requests.get(
            url, stream=True, timeout=(15, 180), headers={"User-Agent": USER_AGENT},
        ) as response:
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "")
            if "text/html" in content_type.lower():
                raise ValueError("GWAS Catalog server returned HTML instead of ZIP data")
            declared = response.headers.get("Content-Length", "")
            if declared.isdigit() and int(declared) > max_download_bytes:
                raise ValueError(
                    f"GWAS Catalog download declares {int(declared):,} bytes, "
                    f"exceeding the {max_download_bytes:,}-byte safety limit"
                )
            with partial.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > max_download_bytes:
                        raise ValueError(
                            "GWAS Catalog download exceeded the configured size limit"
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
                "last_modified": response.headers.get("Last-Modified"),
                "etag": response.headers.get("ETag"),
            }
        # Validate before replacing a previously good cached release.
        downloaded_archive = inspect_archive(partial, max_uncompressed_bytes)
        provenance.update({
            "download_crc_check": downloaded_archive["archive_crc_check"],
            "download_member_name": downloaded_archive["member_name"],
        })
        partial.replace(destination)
        return provenance
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def acquire_archive(args: argparse.Namespace) -> tuple[Path, dict]:
    local = getattr(args, "archive", None)
    if local:
        path = Path(local)
        provenance = {"source_mode": "local", "requested_url": ""}
    else:
        path = Path(args.cache) / "gwas_catalog" / ARCHIVE_FILENAME
        sidecar = path.with_suffix(path.suffix + ".source.json")
        if args.refresh or not path.exists():
            if args.offline:
                raise FileNotFoundError(
                    f"Offline GWAS Catalog cache miss: {path}. Supply --archive or "
                    "run once without --offline."
                )
            provenance = download_archive(
                GWAS_ARCHIVE_URL, path, args.max_download_bytes,
                args.max_uncompressed_bytes,
            )
            provenance["source_mode"] = "downloaded"
            write_json(sidecar, provenance)
        else:
            provenance = {"source_mode": "cache", "requested_url": GWAS_ARCHIVE_URL}
            if sidecar.is_file():
                cached = json.loads(sidecar.read_text(encoding="utf-8"))
                if isinstance(cached, dict):
                    provenance.update(cached)
                    provenance["source_mode"] = "cache"
    inspected = inspect_archive(path, args.max_uncompressed_bytes)
    provenance.update(inspected)
    return path, provenance


def normalize_resolution_rows(rows: Any, rsids: list[str], source: str) -> list[dict]:
    if not isinstance(rows, list):
        raise ValueError(f"dbSNP resolution must be a JSON list: {source}")
    indexed = {str(row.get("input_rsid", "")): row for row in rows
               if isinstance(row, dict)}
    missing = [rsid for rsid in rsids if rsid not in indexed]
    if missing:
        raise ValueError(f"dbSNP resolution is missing panel IDs {missing}: {source}")
    result = []
    for rank, rsid in enumerate(rsids, 1):
        source_row = indexed[rsid]
        current = str(source_row.get("current_rsid", ""))
        if not re.fullmatch(r"rs[0-9]+", current):
            raise ValueError(f"Invalid current_rsid for {rsid} in {source}: {current!r}")
        row = {field: source_row.get(field, "") for field in RESOLUTION_FIELDS}
        row.update({
            "rank": rank,
            "input_rsid": rsid,
            "current_rsid": current,
            "was_merged": current != rsid,
        })
        result.append(row)
    return result


def resolve_dbsnp_ids(args: argparse.Namespace, rsids: list[str], client: Client,
                      warnings: list[dict]) -> tuple[list[dict], str, dict | None]:
    requested = getattr(args, "dbsnp_resolution", None)
    auto_path = DEFAULT_DBSNP_RESOLUTION
    path = Path(requested) if requested else auto_path
    if path.is_file():
        rows = json.loads(path.read_text(encoding="utf-8"))
        return normalize_resolution_rows(rows, rsids, str(path)), str(path.resolve()), None
    if requested:
        raise FileNotFoundError(f"Requested dbSNP resolution file not found: {path}")

    try:
        # Keep the proven resolver in one place until it is promoted to a neutral module.
        from generate_clinvar_dataset import (  # local import avoids startup coupling
            RateLimiter, dbsnp_resolution_row, fetch_dbsnp_summaries,
        )
        records, raw = fetch_dbsnp_summaries(
            client, RateLimiter(0 if args.offline else 0.34), rsids,
        )
        rows = [dbsnp_resolution_row(rsid, rank, records.get(rsid))
                for rank, rsid in enumerate(rsids, 1)]
        if any(not row.get("current_rsid") for row in rows):
            missing = [row["input_rsid"] for row in rows if not row.get("current_rsid")]
            raise ValueError(f"dbSNP returned no current ID for {missing}")
        return rows, "NCBI dbSNP ESummary", raw
    except Exception as exc:
        # Current identifiers are part of the matching contract, not an
        # optional enhancement: this panel includes a known merged rsID. A
        # supplied-ID-only fallback could silently turn a lookup failure into
        # a false "no association" result.
        raise RuntimeError(
            "Unable to resolve every panel rsID to its current dbSNP ID; "
            "supply --dbsnp-resolution or populate the API cache"
        ) from exc


def fetch_metadata(client: Client, enabled: bool, warnings: list[dict]) -> dict:
    if not enabled:
        return {"status": "skipped", "url": GWAS_METADATA_URL}
    try:
        value = client.json(GWAS_METADATA_URL)
        if not isinstance(value, dict):
            raise ValueError("metadata response was not a JSON object")
        return {"status": "available", "url": GWAS_METADATA_URL, **value}
    except Exception as exc:
        warnings.append({
            "severity": "warning", "stage": "gwas_catalog_metadata",
            "message": str(exc),
        })
        return {"status": "unavailable", "url": GWAS_METADATA_URL,
                "error": str(exc)}


def field(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row:
            return row.get(name, "") or ""
    return ""


def relation_row(panel_rsid: str, resolution: dict, rank: int,
                 source_row_number: int, digest: str, row: dict[str, str],
                 reported: list[str], current: list[str], risk: list[str]) -> dict:
    aliases = {panel_rsid, resolution["current_rsid"]}
    reported_matches = [value for value in reported if value in aliases]
    current_matches = [value for value in current if value in aliases]
    if reported_matches and current_matches:
        basis = "reported_and_current_rsid"
    elif reported_matches:
        basis = "reported_rsid"
    else:
        basis = "current_rsid"
    matched = list(dict.fromkeys(reported_matches + current_matches))
    preferred = resolution["current_rsid"]
    matched_rsid = preferred if preferred in matched else matched[0]
    is_multi = len(reported) > 1
    if is_multi:
        scope = "multi_variant_haplotype_or_interaction"
    elif reported:
        scope = "single_variant"
    else:
        scope = "current_id_only"
    return {
        "rank": rank,
        "rsid": panel_rsid,
        "current_rsid": resolution["current_rsid"],
        "was_merged": resolution["current_rsid"] != panel_rsid,
        "source_row_number": source_row_number,
        "source_record_sha256": digest,
        "match_basis": basis,
        "matched_rsid": matched_rsid,
        "matched_identifiers": matched,
        "reported_rsids": reported,
        "current_rsids": current,
        "risk_allele_rsids": risk,
        "association_scope": scope,
        "is_multi_variant": is_multi,
        "strongest_snp_risk_allele": field(row, "STRONGEST SNP-RISK ALLELE"),
        "date_added_to_catalog": field(row, "DATE ADDED TO CATALOG"),
        "pubmed_id": field(row, "PUBMEDID", "PUBMED ID"),
        "first_author": field(row, "FIRST AUTHOR"),
        "publication_date": field(row, "DATE"),
        "journal": field(row, "JOURNAL"),
        "publication_link": field(row, "LINK"),
        "study_title": field(row, "STUDY"),
        "reported_trait": field(row, "DISEASE/TRAIT"),
        "initial_sample_description": field(row, "INITIAL SAMPLE SIZE",
                                             "INITIAL SAMPLE DESCRIPTION"),
        "replication_sample_description": field(row, "REPLICATION SAMPLE SIZE",
                                                 "REPLICATION SAMPLE DESCRIPTION"),
        "region": field(row, "REGION"),
        "chromosome": field(row, "CHR_ID"),
        "position": field(row, "CHR_POS"),
        "reported_genes": field(row, "REPORTED GENE(S)"),
        "mapped_genes": field(row, "MAPPED_GENE", "MAPPED GENE(S)"),
        "context": field(row, "CONTEXT"),
        "intergenic": field(row, "INTERGENIC"),
        "risk_allele_frequency": field(row, "RISK ALLELE FREQUENCY"),
        "p_value": field(row, "P-VALUE"),
        "pvalue_mlog": field(row, "PVALUE_MLOG"),
        "pvalue_text": field(row, "P-VALUE (TEXT)"),
        "or_or_beta": field(row, "OR or BETA"),
        "confidence_interval": field(row, "95% CI (TEXT)"),
        "platform": field(row, "PLATFORM [SNPS PASSING QC]", "PLATFORM (SNPS PASSING QC)"),
        "cnv": field(row, "CNV"),
        "mapped_trait": field(row, "MAPPED_TRAIT"),
        "mapped_trait_uri": field(row, "MAPPED_TRAIT_URI"),
        "study_accession": field(row, "STUDY ACCESSION"),
        "genotyping_technology": field(row, "GENOTYPING TECHNOLOGY"),
        "source_url": f"https://www.ebi.ac.uk/gwas/variants/{resolution['current_rsid']}",
    }


def parse_archive(path: Path, member_name: str, rsids: list[str],
                  resolutions: list[dict]) -> tuple[list[dict], list[dict], list[str], int]:
    resolution_by_input = {row["input_rsid"]: row for row in resolutions}
    aliases: dict[str, list[str]] = {}
    for rsid in rsids:
        for alias in {rsid, resolution_by_input[rsid]["current_rsid"]}:
            aliases.setdefault(alias, []).append(rsid)

    relations: list[dict] = []
    matched_source_rows: list[dict] = []
    header: list[str] = []
    scanned = 0
    with zipfile.ZipFile(path) as archive, archive.open(member_name) as binary:
        with io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = list(reader.fieldnames or [])
            if not header or "SNPS" not in header or "SNP_ID_CURRENT" not in header:
                raise ValueError(
                    "GWAS Catalog association table requires SNPS and "
                    "SNP_ID_CURRENT columns"
                )
            if len(set(header)) != len(header):
                raise ValueError("GWAS Catalog association table has duplicate headers")
            for source_row_number, raw in enumerate(reader, 2):
                scanned += 1
                if None in raw:
                    raise ValueError(
                        f"Malformed GWAS Catalog row {source_row_number}: extra fields"
                    )
                row = {name: raw.get(name, "") or "" for name in header}
                reported = exact_rsid_tokens(row["SNPS"])
                current = current_rsid_tokens(row["SNP_ID_CURRENT"])
                risk = exact_rsid_tokens(field(row, "STRONGEST SNP-RISK ALLELE"))
                matching_inputs: list[str] = []
                for alias in dict.fromkeys(reported + current):
                    for panel_rsid in aliases.get(alias, []):
                        if panel_rsid not in matching_inputs:
                            matching_inputs.append(panel_rsid)
                if not matching_inputs:
                    continue
                digest = source_record_digest(row)
                matched_source_rows.append({
                    "source_row_number": source_row_number,
                    "source_record_sha256": digest,
                    "matched_panel_rsids": sorted(
                        matching_inputs, key=lambda value: rsids.index(value)
                    ),
                    "record": row,
                })
                for panel_rsid in sorted(
                    matching_inputs, key=lambda value: rsids.index(value)
                ):
                    relations.append(relation_row(
                        panel_rsid, resolution_by_input[panel_rsid],
                        rsids.index(panel_rsid) + 1, source_row_number, digest, row,
                        reported, current, risk,
                    ))
    relations.sort(key=lambda row: (row["rank"], row["source_row_number"]))
    return relations, matched_source_rows, header, scanned


def build_statuses(rsids: list[str], resolutions: list[dict],
                   relations: list[dict]) -> list[dict]:
    resolution_by_input = {row["input_rsid"]: row for row in resolutions}
    by_rsid = {rsid: [] for rsid in rsids}
    for relation in relations:
        by_rsid[relation["rsid"]].append(relation)
    statuses = []
    for rank, rsid in enumerate(rsids, 1):
        rows = by_rsid[rsid]
        resolution = resolution_by_input[rsid]
        statuses.append({
            "rank": rank,
            "rsid": rsid,
            "current_rsid": resolution["current_rsid"],
            "was_merged": resolution["current_rsid"] != rsid,
            "dbsnp_status": resolution.get("status", ""),
            "status": "associations_found" if rows else "no_gwas_catalog_association",
            "n_associations": len(rows),
            "n_single_variant_associations": sum(
                row["association_scope"] == "single_variant" for row in rows
            ),
            "n_multi_variant_associations": sum(row["is_multi_variant"] for row in rows),
            "n_matches_by_reported_rsid": sum(
                row["match_basis"] in {"reported_rsid", "reported_and_current_rsid"}
                for row in rows
            ),
            "n_matches_by_current_rsid": sum(
                row["match_basis"] in {"current_rsid", "reported_and_current_rsid"}
                for row in rows
            ),
            "n_unique_studies": len({row["study_accession"] for row in rows
                                      if row["study_accession"]}),
            "n_unique_traits": len({row["mapped_trait"] or row["reported_trait"]
                                     for row in rows
                                     if row["mapped_trait"] or row["reported_trait"]}),
            "error": "",
        })
    return statuses


def sqlite_value(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, bool):
        return int(value)
    return value


def create_table(database: sqlite3.Connection, name: str, fields: list[str],
                 integer_fields: set[str], primary_key: str = "") -> None:
    columns = ", ".join(
        f'"{field_name}" {"INTEGER" if field_name in integer_fields else "TEXT"}'
        for field_name in fields
    )
    suffix = f", PRIMARY KEY ({primary_key})" if primary_key else ""
    database.execute(f'CREATE TABLE "{name}" ({columns}{suffix})')


def insert_rows(database: sqlite3.Connection, name: str, fields: list[str],
                rows: Iterable[dict]) -> None:
    columns = ", ".join(f'"{value}"' for value in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{name}" ({columns}) VALUES ({placeholders})',
        [tuple(sqlite_value(row.get(value, "")) for value in fields) for row in rows],
    )


def write_sqlite(path: Path, resolutions: list[dict], relations: list[dict],
                 raw_rows: list[dict], statuses: list[dict], manifest: dict) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    raw_fields = ["source_row_number", "source_record_sha256",
                  "matched_panel_rsids", "record"]
    try:
        with database:
            create_table(database, "dbsnp_resolution", RESOLUTION_FIELDS,
                         {"rank", "was_merged"}, '"input_rsid"')
            create_table(database, "associations", ASSOCIATION_FIELDS,
                         {"rank", "was_merged", "source_row_number", "is_multi_variant"},
                         '"rsid", "source_row_number"')
            create_table(database, "matched_source_records", raw_fields,
                         {"source_row_number"}, '"source_row_number"')
            status_integers = {"rank", "was_merged"} | {
                name for name in STATUS_FIELDS if name.startswith("n_")
            }
            create_table(database, "panel_status", STATUS_FIELDS, status_integers,
                         '"rsid"')
            database.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)")
            insert_rows(database, "dbsnp_resolution", RESOLUTION_FIELDS, resolutions)
            insert_rows(database, "associations", ASSOCIATION_FIELDS, relations)
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
    text = f"""# GWAS Catalog dataset for the QIAGEN 20-rsID panel

This folder is a reproducible snapshot of literature-curated top associations
from the NHGRI-EBI GWAS Catalog for the exact 20 public test rsIDs. It contains
no patient genotypes.

- Metadata-reported current release: {manifest.get('gwas_catalog_data_release', 'not available')}
- Archive SHA-256: {manifest['archive_provenance']['archive_sha256']}
- Input rsIDs: {manifest['n_input_rsids']}
- rsIDs with one or more associations: {manifest['n_rsids_with_associations']}
- Normalized rsID-association relationships: {manifest['n_association_relations']}
- Losslessly retained source rows: {manifest['n_matched_source_rows']}
- Multi-variant/haplotype-or-interaction relationships: {manifest['n_multi_variant_relations']}

## Files

`rsid_associations.tsv` is the analysis-friendly table with one relationship
between a panel rsID and a source association row. `matched_source_records.json`
retains every matching official row and all original strings without numeric
coercion; its TSV copy adds only provenance and match columns. A source row can
produce more than one relationship when it contains multiple panel rsIDs.

`panel_status.tsv` always contains all 20 inputs in their original order.
`dbsnp_resolution.tsv` records current identifiers used during matching.
`gwas_catalog_metadata.json`, `archive_provenance.json`, and `run_manifest.json`
capture metadata and file provenance. The metadata release is recorded
separately and is not asserted to identify a user-supplied or previously
cached archive; the archive hash is authoritative. `gwas_catalog.sqlite` mirrors the main
tables. `errors.json` records non-fatal lookup warnings.

`association_scope=multi_variant_haplotype_or_interaction` means the official
`SNPS` field contains multiple rsIDs. Such a set-level association must not be
assigned to one component SNP alone. Reported risk/effect alleles are retained
as published; they are not assumed to be forward-strand GRCh38 alternate
alleles and are not harmonized here.

GWAS association is not proof of causality, a diagnosis, or an individual risk
estimate. This archive contains curated top associations, not every result from
full genome-wide summary statistics.

Official sources:
- https://www.ebi.ac.uk/gwas/docs/programmatic-access/rest-api/
- https://www.ebi.ac.uk/gwas/docs/file-downloads/
- https://www.ebi.ac.uk/gwas/docs/fileheaders/
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def run(args: argparse.Namespace, client: Client | None = None) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    if args.max_download_bytes < 1 or args.max_uncompressed_bytes < 1:
        raise ValueError("Archive size limits must be positive")
    rsids = select_panel_rsids(args.input, require_exact=True)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # Never let a failed refresh leave an older generated database or raw
    # resolver response looking current. Only files owned by this generator
    # are removed; unrelated files in the output directory are preserved.
    for stale_name in ("gwas_catalog.sqlite", "gwas_catalog.sqlite.part"):
        (out / stale_name).unlink(missing_ok=True)
    (out / "raw" / "dbsnp_esummary.json").unlink(missing_ok=True)
    warnings: list[dict] = []
    client = client or Client(args.cache, args.offline, args.refresh)

    archive_path, archive_provenance = acquire_archive(args)
    metadata = fetch_metadata(client, not args.no_metadata, warnings)
    resolutions, resolution_source, dbsnp_raw = resolve_dbsnp_ids(
        args, rsids, client, warnings,
    )
    if dbsnp_raw is not None:
        write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)

    relations, raw_rows, source_header, scanned = parse_archive(
        archive_path, archive_provenance["member_name"], rsids, resolutions,
    )
    statuses = build_statuses(rsids, resolutions, relations)

    input_path = Path(args.input)
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source_database": "NHGRI-EBI GWAS Catalog",
        "source_archive_url": archive_provenance.get("requested_url", ""),
        "gwas_catalog_data_release": metadata.get("data_release_date", ""),
        "gwas_catalog_release_scope": (
            "REST metadata snapshot; not cryptographically bound to the archive"
        ),
        "gwas_catalog_genome_build": metadata.get("gene_build", ""),
        "gwas_catalog_dbsnp_build": metadata.get("dbsnp_build", ""),
        "dbsnp_resolution_source": resolution_source,
        "n_input_rsids": len(rsids),
        "n_merged_input_rsids": sum(bool(row["was_merged"]) for row in resolutions),
        "n_source_rows_scanned": scanned,
        "n_matched_source_rows": len(raw_rows),
        "n_association_relations": len(relations),
        "n_multi_variant_relations": sum(row["is_multi_variant"] for row in relations),
        "n_rsids_with_associations": sum(row["n_associations"] > 0 for row in statuses),
        "n_warnings": len(warnings),
        "source_header": source_header,
        "archive_provenance": archive_provenance,
        "http_provenance": list(client.manifest.values()),
        "notes": [
            "Exact rsID tokens are matched against supplied and current dbSNP identifiers",
            "All original source values remain strings in matched_source_records",
            "Multi-rsID rows are retained and explicitly labeled",
            "Risk/effect alleles are reported values and are not assumed to be genomic ALT",
            "A GWAS association is not evidence that the matched variant is causal",
            "The REST metadata call is optional; the official bulk archive is the data source",
        ],
    }

    raw_tsv_rows = [{
        "source_row_number": item["source_row_number"],
        "source_record_sha256": item["source_record_sha256"],
        "matched_panel_rsids": item["matched_panel_rsids"],
        **item["record"],
    } for item in raw_rows]
    raw_tsv_fields = ["source_row_number", "source_record_sha256",
                      "matched_panel_rsids", *source_header]

    write_json(out / "archive_provenance.json", archive_provenance)
    write_json(out / "gwas_catalog_metadata.json", metadata)
    write_json(out / "dbsnp_resolution.json", resolutions)
    write_tsv(out / "dbsnp_resolution.tsv", resolutions, RESOLUTION_FIELDS)
    write_json(out / "matched_source_records.json", raw_rows)
    write_tsv(out / "matched_source_records.tsv", raw_tsv_rows, raw_tsv_fields)
    write_json(out / "rsid_associations.json", relations)
    write_tsv(out / "rsid_associations.tsv", relations, ASSOCIATION_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "errors.json", warnings)
    write_json(out / "run_manifest.json", manifest)
    write_sqlite(out / "gwas_catalog.sqlite", resolutions, relations, raw_rows,
                 statuses, manifest)
    write_readme(out, manifest)
    print(
        f"GWAS Catalog dataset complete: {len(rsids)} input rsIDs, "
        f"{len(relations)} association relationships from {len(raw_rows)} source rows. "
        f"See {out}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="Exact 20-rsID QIAGEN TSV in original order")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument("--archive", help="Local official association ZIP (no download)")
    parser.add_argument("--cache", default=".cache", help="Download/API cache directory")
    parser.add_argument("--offline", action="store_true",
                        help="Use --archive or the cached ZIP and cached API data only")
    parser.add_argument("--refresh", action="store_true",
                        help="Refresh the official ZIP and API cache")
    parser.add_argument(
        "--dbsnp-resolution",
        help=("Existing dbsnp_resolution.json; by default reuse the sibling "
              "ClinVar dataset, then try NCBI ESummary"),
    )
    parser.add_argument("--no-metadata", action="store_true",
                        help="Skip the optional GWAS Catalog v2 metadata request")
    parser.add_argument("--max-download-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--max-uncompressed-bytes", type=int, default=4_000_000_000)
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
