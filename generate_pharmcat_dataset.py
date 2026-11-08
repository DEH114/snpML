"""Build a PharmCAT reference-position dataset for an rsID panel.

This generator matches supplied and current dbSNP identifiers against the
official PharmCAT GRCh38 reference-position VCF.  It reports locus coverage
and the source ``PX``/``POI`` annotations only.  It does not run PharmCAT,
interpret the template genotype column, infer star alleles or diplotypes, or
produce prescribing recommendations.  Those operations require a validated
patient/sample VCF and are outside this rsID-only workflow.
"""
from __future__ import annotations

import argparse
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
from generate_gwas_catalog_dataset import RESOLUTION_FIELDS, resolve_dbsnp_ids


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "pharmcat_dataset"
PHARMCAT_VERSION = "3.4.0"
PHARMCAT_TAG = f"v{PHARMCAT_VERSION}"
POSITIONS_FILENAME = "pharmcat_positions.vcf"
POSITIONS_URL = (
    "https://raw.githubusercontent.com/PharmGKB/PharmCAT/"
    f"{PHARMCAT_TAG}/{POSITIONS_FILENAME}"
)
PHARMCAT_REPOSITORY_URL = "https://github.com/PharmGKB/PharmCAT"
USER_AGENT = "BioHack-SNP-Toolkit/1.0 (research prototype)"
DEFAULT_MAX_DOWNLOAD_BYTES = 10 * 1024 * 1024

FIXED_VCF_FIELDS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
EXPECTED_VCF_FIELDS = [*FIXED_VCF_FIELDS, "FORMAT", "PharmCAT"]
RSID_PATTERN = re.compile(r"rs[1-9][0-9]*")

MATCH_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "source_line_number",
    "source_record_sha256", "match_basis", "matched_rsid", "chromosome",
    "position", "source_id", "reference", "alternate_alleles", "quality",
    "filter", "info_raw", "pharmacogenes", "position_of_interest",
    "allele_definition_position", "source_url",
]

STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "dbsnp_status", "status",
    "n_position_matches", "n_allele_definition_positions",
    "n_positions_of_interest", "pharmacogenes", "matched_source_rsids", "error",
]

OWNED_OUTPUTS = (
    "source_provenance.json", "dbsnp_resolution.json", "dbsnp_resolution.tsv",
    "matched_source_records.json", "matched_source_records.tsv",
    "pharmcat_position_matches.json", "pharmcat_position_matches.tsv",
    "panel_status.json", "panel_status.tsv", "errors.json", "run_manifest.json",
    "pharmcat.sqlite", "pharmcat.sqlite.part", "README.md",
)


def normalize_release(value: str) -> tuple[str, str]:
    """Return the numeric release and canonical Git tag for an official release."""
    supplied = str(value).strip()
    version = supplied[1:] if supplied.startswith("v") else supplied
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise ValueError(
            "--release must be an official semantic version such as 3.4.0 or v3.4.0"
        )
    return version, f"v{version}"


def positions_url_for_tag(tag: str) -> str:
    return (
        "https://raw.githubusercontent.com/PharmGKB/PharmCAT/"
        f"{tag}/{POSITIONS_FILENAME}"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_line_ending(text: str) -> tuple[str, str]:
    if text.endswith("\r\n"):
        return text[:-2], "\r\n"
    if text.endswith(("\n", "\r")):
        return text[:-1], text[-1]
    return text, ""


def _read_utf8_lines(path: Path) -> list[tuple[int, bytes, str, str]]:
    try:
        raw_lines = path.read_bytes().splitlines(keepends=True)
    except OSError as exc:
        raise ValueError(f"Cannot read PharmCAT positions VCF {path}: {exc}") from exc
    result = []
    for line_number, raw_line in enumerate(raw_lines, 1):
        try:
            decoded = raw_line.decode("utf-8-sig" if line_number == 1 else "utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"PharmCAT positions VCF is not valid UTF-8 at line {line_number}"
            ) from exc
        line, ending = _split_line_ending(decoded)
        result.append((line_number, raw_line, line, ending))
    return result


def inspect_positions_vcf(path: Path, max_bytes: int) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"PharmCAT positions VCF not found: {path}")
    size = path.stat().st_size
    if size < 1:
        raise ValueError("PharmCAT positions VCF is empty")
    if size > max_bytes:
        raise ValueError(
            f"PharmCAT positions VCF is {size:,} bytes, exceeding the "
            f"{max_bytes:,}-byte safety limit"
        )
    with path.open("rb") as handle:
        prefix = handle.read(4096).lstrip().lower()
    if prefix.startswith((b"<!doctype html", b"<html", b"<?xml")):
        raise ValueError("PharmCAT server returned HTML/XML instead of a VCF")

    lines = _read_utf8_lines(path)
    metadata: list[str] = []
    header: list[str] | None = None
    header_line_number = 0
    for line_number, _raw, line, _ending in lines:
        if not line:
            raise ValueError(f"PharmCAT positions VCF has a blank line at {line_number}")
        if line.startswith("##"):
            if header is not None:
                raise ValueError("PharmCAT positions VCF metadata appears after its header")
            metadata.append(line)
            continue
        if line.startswith("#"):
            if header is not None:
                raise ValueError("PharmCAT positions VCF has duplicate column headers")
            header = line.split("\t")
            header_line_number = line_number
            continue
        if header is None:
            raise ValueError("PharmCAT positions VCF has data before its column header")

    if "##fileformat=VCFv4.2" not in metadata:
        raise ValueError("PharmCAT positions VCF must declare ##fileformat=VCFv4.2")
    source_values = [line for line in metadata if line.startswith("##source=")]
    if not source_values or not any("pharmcat" in line.casefold() for line in source_values):
        raise ValueError("PharmCAT positions VCF has no PharmCAT source declaration")
    contigs = [line for line in metadata if line.startswith("##contig=<")]
    if not contigs:
        raise ValueError("PharmCAT positions VCF has no contig metadata")
    invalid_assemblies = [line for line in contigs if "assembly=GRCh38" not in line]
    if invalid_assemblies:
        raise ValueError("PharmCAT positions VCF must use GRCh38 contig metadata")
    if not any(line.startswith("##INFO=<ID=PX,") for line in metadata):
        raise ValueError("PharmCAT positions VCF is missing its PX INFO declaration")
    if not any(line.startswith("##INFO=<ID=POI,") for line in metadata):
        raise ValueError("PharmCAT positions VCF is missing its POI INFO declaration")
    if header != EXPECTED_VCF_FIELDS:
        raise ValueError(
            "PharmCAT positions VCF has an unexpected column header: "
            f"expected={EXPECTED_VCF_FIELDS}, observed={header}"
        )
    return {
        "positions_vcf_path": str(path.resolve()),
        "positions_vcf_bytes": size,
        "positions_vcf_sha256": sha256_file(path),
        "vcf_header": header,
        "vcf_header_line_number": header_line_number,
        "vcf_metadata": metadata,
        "genome_assembly": "GRCh38",
    }


def download_positions_vcf(url: str, destination: Path, max_bytes: int) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    partial.unlink(missing_ok=True)
    digest = hashlib.sha256()
    total = 0
    try:
        with requests.get(
            url, stream=True, timeout=(15, 180),
            headers={"User-Agent": USER_AGENT, "Accept": "text/plain"},
        ) as response:
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "")
            if "text/html" in content_type.casefold():
                raise ValueError("PharmCAT server returned HTML instead of a VCF")
            declared = response.headers.get("Content-Length", "")
            if declared.isdigit() and int(declared) > max_bytes:
                raise ValueError(
                    f"PharmCAT download declares {int(declared):,} bytes, exceeding "
                    f"the {max_bytes:,}-byte safety limit"
                )
            with partial.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > max_bytes:
                        raise ValueError("PharmCAT download exceeded the configured size limit")
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
        inspect_positions_vcf(partial, max_bytes)
        partial.replace(destination)
        return provenance
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def acquire_positions_vcf(args: argparse.Namespace) -> tuple[Path, dict[str, Any]]:
    version, tag = normalize_release(
        getattr(args, "release", PHARMCAT_VERSION) or PHARMCAT_VERSION
    )
    positions_url = positions_url_for_tag(tag)
    local = getattr(args, "positions_vcf", None)
    if local:
        path = Path(local)
        provenance: dict[str, Any] = {"source_mode": "local", "requested_url": ""}
    else:
        path = Path(args.cache) / "pharmcat" / version / POSITIONS_FILENAME
        sidecar = path.with_suffix(path.suffix + ".source.json")
        if args.refresh or not path.exists():
            if args.offline:
                raise FileNotFoundError(
                    f"Offline PharmCAT cache miss: {path}. Supply --positions-vcf "
                    "or run once without --offline."
                )
            provenance = download_positions_vcf(
                positions_url, path, args.max_download_bytes,
            )
            provenance["source_mode"] = "downloaded"
            write_json(sidecar, provenance)
        else:
            provenance = {"source_mode": "cache", "requested_url": positions_url}
            if sidecar.is_file():
                try:
                    cached = json.loads(sidecar.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError) as exc:
                    raise ValueError(f"Invalid PharmCAT cache provenance: {sidecar}") from exc
                if not isinstance(cached, dict):
                    raise ValueError(f"PharmCAT cache provenance must be an object: {sidecar}")
                provenance.update(cached)
                provenance["source_mode"] = "cache"
    provenance.update(inspect_positions_vcf(path, args.max_download_bytes))
    provenance.update({
        "pharmcat_version": version,
        "pharmcat_tag": tag,
        "official_positions_url": positions_url,
    })
    return path, provenance


def parse_info(value: str, line_number: int) -> dict[str, str | bool]:
    if value in {"", "."}:
        return {}
    result: dict[str, str | bool] = {}
    for item in value.split(";"):
        if not item:
            raise ValueError(f"Malformed empty INFO item at VCF line {line_number}")
        if "=" in item:
            key, item_value = item.split("=", 1)
            if not key or item_value == "":
                raise ValueError(f"Malformed INFO item {item!r} at VCF line {line_number}")
            value_to_store: str | bool = item_value
        else:
            key, value_to_store = item, True
        if key in result:
            raise ValueError(f"Duplicate INFO key {key!r} at VCF line {line_number}")
        result[key] = value_to_store
    return result


def source_ids(value: str, line_number: int) -> list[str]:
    if value == ".":
        return []
    values = value.split(";")
    if any(not RSID_PATTERN.fullmatch(item) for item in values):
        raise ValueError(
            f"PharmCAT VCF line {line_number} contains a non-canonical ID: {value!r}"
        )
    if len(values) != len(set(values)):
        raise ValueError(f"PharmCAT VCF line {line_number} contains duplicate IDs")
    return values


def normalized_match(
    panel_rsid: str, resolution: dict, rank: int, line_number: int,
    digest: str, row: dict[str, str], matched_rsid: str,
    info: dict[str, str | bool], source_url: str,
) -> dict[str, Any]:
    current = str(resolution["current_rsid"])
    match_basis = (
        "current_rsid" if matched_rsid == current and current != panel_rsid
        else "input_rsid"
    )
    px = info.get("PX", "")
    pharmacogenes = list(dict.fromkeys(
        gene.strip() for gene in str(px).split(",") if gene.strip()
    )) if isinstance(px, str) else []
    alts = row["ALT"].split(",")
    return {
        "rank": rank,
        "rsid": panel_rsid,
        "current_rsid": current,
        "was_merged": current != panel_rsid,
        "source_line_number": line_number,
        "source_record_sha256": digest,
        "match_basis": match_basis,
        "matched_rsid": matched_rsid,
        "chromosome": row["#CHROM"],
        "position": int(row["POS"]),
        "source_id": row["ID"],
        "reference": row["REF"],
        "alternate_alleles": alts,
        "quality": row["QUAL"],
        "filter": row["FILTER"],
        "info_raw": row["INFO"],
        "pharmacogenes": pharmacogenes,
        "position_of_interest": "POI" in info,
        "allele_definition_position": bool(pharmacogenes),
        "source_url": source_url,
    }


def parse_positions_vcf(
    path: Path, rsids: list[str], resolutions: list[dict],
    source_url: str = POSITIONS_URL,
) -> tuple[list[dict], list[dict], list[str], int]:
    inspected = inspect_positions_vcf(path, max(path.stat().st_size, 1))
    header = list(inspected["vcf_header"])
    header_line = int(inspected["vcf_header_line_number"])
    resolution_by_input = {str(row["input_rsid"]): row for row in resolutions}
    if list(resolution_by_input) != rsids:
        raise ValueError("dbSNP resolution does not cover the exact input panel in order")
    rank_by_rsid = {rsid: rank for rank, rsid in enumerate(rsids, 1)}
    aliases: dict[str, list[str]] = {}
    for rsid in rsids:
        current = str(resolution_by_input[rsid]["current_rsid"])
        for alias in dict.fromkeys((rsid, current)):
            aliases.setdefault(alias, []).append(rsid)

    matches: list[dict] = []
    matched_records: list[dict] = []
    scanned = 0
    for line_number, raw_line, source_line, line_ending in _read_utf8_lines(path):
        if line_number <= header_line:
            continue
        if source_line.startswith("#"):
            raise ValueError(f"Unexpected header/comment at VCF line {line_number}")
        scanned += 1
        values = source_line.split("\t")
        if len(values) != len(header):
            raise ValueError(
                f"Malformed PharmCAT VCF line {line_number}: expected "
                f"{len(header)} fields, observed {len(values)}"
            )
        row = dict(zip(header, values))
        if not re.fullmatch(r"[1-9][0-9]*", row["POS"]):
            raise ValueError(f"Invalid VCF position at line {line_number}: {row['POS']!r}")
        if not row["#CHROM"].startswith("chr"):
            raise ValueError(f"Invalid GRCh38 chromosome at VCF line {line_number}")
        if not row["REF"] or row["REF"] == ".":
            raise ValueError(f"Missing REF at VCF line {line_number}")
        alts = row["ALT"].split(",")
        if any(not alt or alt == "." for alt in alts):
            raise ValueError(f"Missing ALT at VCF line {line_number}")
        if row["FORMAT"] != "GT":
            raise ValueError(
                f"PharmCAT template FORMAT must be GT at VCF line {line_number}"
            )
        ids = source_ids(row["ID"], line_number)
        info = parse_info(row["INFO"], line_number)
        if "PX" in info and not isinstance(info["PX"], str):
            raise ValueError(f"PX must have a value at VCF line {line_number}")
        if "POI" in info and info["POI"] is not True:
            raise ValueError(f"POI must be a flag at VCF line {line_number}")
        matching: dict[str, str] = {}
        for source_rsid in ids:
            for panel_rsid in aliases.get(source_rsid, []):
                matching.setdefault(panel_rsid, source_rsid)
        if not matching:
            continue
        ordered_inputs = sorted(matching, key=rank_by_rsid.__getitem__)
        digest = hashlib.sha256(raw_line).hexdigest()
        matched_records.append({
            "source_line_number": line_number,
            "source_record_sha256": digest,
            "matched_panel_rsids": ordered_inputs,
            "source_line": source_line,
            "line_ending": line_ending,
            "record": row,
        })
        for panel_rsid in ordered_inputs:
            matches.append(normalized_match(
                panel_rsid, resolution_by_input[panel_rsid],
                rank_by_rsid[panel_rsid], line_number, digest, row,
                matching[panel_rsid], info, source_url,
            ))
    matches.sort(key=lambda row: (row["rank"], row["source_line_number"]))
    return matches, matched_records, header, scanned


def build_statuses(
    rsids: list[str], resolutions: list[dict], matches: list[dict],
) -> list[dict]:
    resolution_by_input = {str(row["input_rsid"]): row for row in resolutions}
    by_rsid = {rsid: [] for rsid in rsids}
    for match in matches:
        by_rsid[str(match["rsid"])].append(match)
    statuses = []
    for rank, rsid in enumerate(rsids, 1):
        rows = by_rsid[rsid]
        resolution = resolution_by_input[rsid]
        statuses.append({
            "rank": rank,
            "rsid": rsid,
            "current_rsid": resolution["current_rsid"],
            "was_merged": bool(resolution["was_merged"]),
            "dbsnp_status": resolution.get("status", ""),
            "status": "pharmcat_positions_found" if rows else "no_pharmcat_position",
            "n_position_matches": len(rows),
            "n_allele_definition_positions": sum(
                bool(row["allele_definition_position"]) for row in rows
            ),
            "n_positions_of_interest": sum(bool(row["position_of_interest"]) for row in rows),
            "pharmacogenes": list(dict.fromkeys(
                gene for row in rows for gene in row["pharmacogenes"]
            )),
            "matched_source_rsids": list(dict.fromkeys(
                row["matched_rsid"] for row in rows
            )),
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
        f'"{field}" {"INTEGER" if field in integer_fields else "TEXT"}'
        for field in fields
    )
    suffix = f", PRIMARY KEY ({primary_key})" if primary_key else ""
    database.execute(f'CREATE TABLE "{name}" ({columns}{suffix})')


def insert_rows(database: sqlite3.Connection, name: str, fields: list[str],
                rows: Iterable[dict]) -> None:
    columns = ", ".join(f'"{field}"' for field in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{name}" ({columns}) VALUES ({placeholders})',
        [tuple(sqlite_value(row.get(field, "")) for field in fields) for row in rows],
    )


def write_sqlite(path: Path, resolutions: list[dict], matches: list[dict],
                 source_rows: list[dict], statuses: list[dict], manifest: dict) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    raw_fields = [
        "source_line_number", "source_record_sha256", "matched_panel_rsids",
        "source_line", "line_ending", "record",
    ]
    try:
        with database:
            create_table(database, "dbsnp_resolution", RESOLUTION_FIELDS,
                         {"rank", "was_merged"}, '"input_rsid"')
            create_table(database, "pharmcat_position_matches", MATCH_FIELDS,
                         {"rank", "was_merged", "source_line_number", "position",
                          "position_of_interest", "allele_definition_position"},
                         '"rsid", "source_line_number"')
            create_table(database, "matched_source_records", raw_fields,
                         {"source_line_number"}, '"source_line_number"')
            status_integers = {"rank", "was_merged"} | {
                field for field in STATUS_FIELDS if field.startswith("n_")
            }
            create_table(database, "panel_status", STATUS_FIELDS, status_integers,
                         '"rsid"')
            database.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)"
            )
            insert_rows(database, "dbsnp_resolution", RESOLUTION_FIELDS, resolutions)
            insert_rows(database, "pharmcat_position_matches", MATCH_FIELDS, matches)
            insert_rows(database, "matched_source_records", raw_fields, source_rows)
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
    text = f"""# PharmCAT reference-position dataset

This folder records which loci in an ordered {manifest['n_input_rsids']}-rsID
panel occur in the official PharmCAT {manifest['pharmcat_version']} GRCh38
reference-position VCF.

- Source VCF SHA-256: {manifest['source_provenance']['positions_vcf_sha256']}
- Input rsIDs: {manifest['n_input_rsids']}
- rsIDs with one or more reference-position matches: {manifest['n_rsids_with_position_matches']}
- Normalized position relationships: {manifest['n_position_match_relations']}
- Losslessly retained source rows: {manifest['n_matched_source_rows']}

`pharmcat_position_matches.tsv` contains exact supplied/current-rsID matches and
the source `PX` pharmacogene and `POI` flag. `matched_source_records.json`
preserves the complete original UTF-8 VCF source line, its exact line ending,
all fields, and a SHA-256 of the original bytes. `panel_status.tsv` contains
every input in original order. `pharmcat.sqlite` mirrors these normalized
tables, and the provenance and run manifests identify the pinned source.

## Scope limitation

This is a **locus-only reference-position index**, not a PharmCAT clinical
report. The `PharmCAT` sample column in the official reference VCF is a template
and is retained only for source fidelity; this generator never interprets its
`GT` value. An rsID list provides no patient genotype, phase, copy number, or
quality evidence. Consequently this dataset does not infer star alleles,
diplotypes, metabolizer phenotypes, drug response, or prescribing
recommendations. A missing position is not evidence of no pharmacogenomic
importance.

Official sources:
- https://github.com/PharmGKB/PharmCAT/releases/tag/{manifest['pharmcat_tag']}
- {PHARMCAT_REPOSITORY_URL}
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
    for filename in OWNED_OUTPUTS:
        (out / filename).unlink(missing_ok=True)
    (out / "raw" / "dbsnp_esummary.json").unlink(missing_ok=True)

    warnings: list[dict] = []
    client = client or Client(args.cache, args.offline, args.refresh)
    positions_path, source_provenance = acquire_positions_vcf(args)
    resolutions, resolution_source, dbsnp_raw = resolve_dbsnp_ids(
        args, rsids, client, warnings,
    )
    if dbsnp_raw is not None:
        write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)

    matches, source_rows, source_header, scanned = parse_positions_vcf(
        positions_path, rsids, resolutions,
        str(source_provenance["official_positions_url"]),
    )
    statuses = build_statuses(rsids, resolutions, matches)
    input_path = Path(args.input)
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "source_database": "PharmCAT reference positions",
        "pharmcat_version": source_provenance["pharmcat_version"],
        "pharmcat_tag": source_provenance["pharmcat_tag"],
        "source_positions_url": source_provenance["official_positions_url"],
        "genome_assembly": "GRCh38",
        "dbsnp_resolution_source": resolution_source,
        "n_input_rsids": len(rsids),
        "n_merged_input_rsids": sum(bool(row["was_merged"]) for row in resolutions),
        "n_source_rows_scanned": scanned,
        "n_matched_source_rows": len(source_rows),
        "n_position_match_relations": len(matches),
        "n_rsids_with_position_matches": sum(
            row["n_position_matches"] > 0 for row in statuses
        ),
        "n_allele_definition_position_relations": sum(
            bool(row["allele_definition_position"]) for row in matches
        ),
        "n_position_of_interest_relations": sum(
            bool(row["position_of_interest"]) for row in matches
        ),
        "n_warnings": len(warnings),
        "source_header": source_header,
        "source_provenance": source_provenance,
        "http_provenance": list(client.manifest.values()),
        "notes": [
            "Exact supplied and current dbSNP identifiers are matched",
            "PX and POI are retained from the pinned PharmCAT reference-position VCF",
            "Matched VCF source lines and line endings are preserved losslessly",
            "The PharmCAT template genotype column is never interpreted",
            "An rsID list contains no patient genotype, phase, copy number, or quality",
            "No star allele, diplotype, phenotype, or prescribing recommendation is inferred",
        ],
    }

    raw_tsv_rows = [{
        "source_line_number": item["source_line_number"],
        "source_record_sha256": item["source_record_sha256"],
        "matched_panel_rsids": item["matched_panel_rsids"],
        "source_line": item["source_line"],
        "line_ending": item["line_ending"].encode("unicode_escape").decode("ascii"),
        **item["record"],
    } for item in source_rows]
    raw_tsv_fields = [
        "source_line_number", "source_record_sha256", "matched_panel_rsids",
        "source_line", "line_ending", *source_header,
    ]

    write_json(out / "source_provenance.json", source_provenance)
    write_json(out / "dbsnp_resolution.json", resolutions)
    write_tsv(out / "dbsnp_resolution.tsv", resolutions, RESOLUTION_FIELDS)
    write_json(out / "matched_source_records.json", source_rows)
    write_tsv(out / "matched_source_records.tsv", raw_tsv_rows, raw_tsv_fields)
    write_json(out / "pharmcat_position_matches.json", matches)
    write_tsv(out / "pharmcat_position_matches.tsv", matches, MATCH_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "errors.json", warnings)
    write_json(out / "run_manifest.json", manifest)
    write_sqlite(
        out / "pharmcat.sqlite", resolutions, matches, source_rows, statuses, manifest,
    )
    write_readme(out, manifest)
    print(
        f"PharmCAT reference-position dataset complete: {len(rsids)} input rsIDs, "
        f"{len(matches)} position relationships from {len(source_rows)} source rows. "
        f"See {out}"
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
        "--release", default=PHARMCAT_VERSION,
        help=("Official PharmCAT release tag/version used for the positions VCF "
              f"(default: {PHARMCAT_VERSION})"),
    )
    parser.add_argument(
        "--positions-vcf", help="Local official PharmCAT positions VCF (no download)",
    )
    parser.add_argument("--cache", default=".cache", help="Download/API cache directory")
    parser.add_argument(
        "--offline", action="store_true",
        help="Use --positions-vcf or the pinned cached VCF and cached API data only",
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Refresh the pinned official VCF and API cache",
    )
    parser.add_argument(
        "--dbsnp-resolution",
        help=("Existing dbsnp_resolution.json; by default reuse the sibling "
              "ClinVar dataset, then try NCBI ESummary"),
    )
    parser.add_argument(
        "--max-download-bytes", type=int, default=DEFAULT_MAX_DOWNLOAD_BYTES,
        help="Maximum accepted downloaded, cached, or local VCF size",
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
