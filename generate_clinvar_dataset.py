"""Build a ClinVar database for the exact 20 QIAGEN challenge rsIDs.

The generator first resolves current/merged identifiers with NCBI dbSNP,
follows dbSNP's explicit links to ClinVar, and then downloads ClinVar ESummary
records. Exact current dbSNP cross-references are verified in every linked
ClinVar record. Direct variants are kept separate from haplotype/genotype
records so a set-level interpretation is not assigned to one component SNP.
Outputs include TSV, JSON, raw API responses, and a SQLite database.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import Client, utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "clinvar_dataset"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "BioHack_SNP_Toolkit"
SUMMARY_CHUNK_SIZE = 100

GRCH38_ACCESSION_TO_CHROM = {
    **{f"NC_{number:06d}": str(number) for number in range(1, 23)},
    "NC_000023": "X",
    "NC_000024": "Y",
    "NC_012920": "MT",
}

DBSNP_FIELDS = [
    "rank", "input_rsid", "current_rsid", "was_merged", "dbsnp_uid",
    "snp_class", "spdi", "chromosome", "chrpos_zero_based", "merge_text",
    "status", "error",
]
RECORD_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "clinvar_link_rank",
    "clinvar_uid", "accession", "accession_version", "record_scope",
    "exact_current_dbsnp_xref", "matching_allele_roles", "object_type", "title", "variant_keys",
    "matching_components", "component_rsids", "genes",
    "molecular_consequences", "protein_change",
    "germline_classification", "germline_review_status",
    "germline_last_evaluated", "germline_traits",
    "clinical_impact_classification", "clinical_impact_review_status",
    "clinical_impact_last_evaluated", "clinical_impact_traits",
    "oncogenicity_classification", "oncogenicity_review_status",
    "oncogenicity_last_evaluated", "oncogenicity_traits",
    "n_scv", "scv_accessions", "n_rcv", "rcv_accessions", "source_url",
]
LINK_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "clinvar_link_rank",
    "clinvar_uid", "accession", "accession_version", "record_scope",
    "allele_role", "matching_variant_keys", "matching_component_count",
    "exact_current_dbsnp_xref", "source_url",
]
STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "dbsnp_status",
    "dbsnp_snp_class", "dbsnp_spdi", "status", "n_clinvar_links",
    "n_linked_summaries", "n_exact_records", "n_direct_records",
    "n_related_set_records", "n_rejected_linked_records",
    "direct_accessions", "direct_germline_classifications",
    "direct_review_statuses", "error",
]
TRAIT_FIELDS = [
    "rsid", "current_rsid", "clinvar_uid", "record_scope",
    "classification_type", "classification_description", "review_status",
    "trait_name", "trait_xrefs",
]


class RateLimiter:
    """Keep NCBI requests at or below the unauthenticated three/second limit."""

    def __init__(self, interval: float = 0.34):
        self.interval = interval
        self.last_request = 0.0

    def wait(self) -> None:
        delay = self.interval - (time.monotonic() - self.last_request)
        if delay > 0:
            time.sleep(delay)
        self.last_request = time.monotonic()


def request_json(client: Client, limiter: RateLimiter, endpoint: str,
                 params: dict[str, Any]) -> dict:
    limiter.wait()
    data = client.json(f"{EUTILS}/{endpoint}", params)
    if not isinstance(data, dict):
        raise ValueError(f"NCBI {endpoint} returned a non-object JSON response")
    return data


def chunks(values: list[str], size: int):
    for start in range(0, len(values), size):
        yield values[start:start + size]


def fetch_dbsnp_summaries(client: Client, limiter: RateLimiter,
                          rsids: list[str]) -> tuple[dict[str, dict], dict]:
    """Fetch dbSNP summaries keyed by the supplied (possibly merged) rsID."""
    records: dict[str, dict] = {}
    responses = []
    numeric_ids = [rsid.removeprefix("rs") for rsid in rsids]
    for group in chunks(numeric_ids, SUMMARY_CHUNK_SIZE):
        response = request_json(client, limiter, "esummary.fcgi", {
            "db": "snp", "id": ",".join(group), "retmode": "json",
            "version": "2.0", "tool": TOOL,
        })
        responses.append(response)
        result = response.get("result", {})
        for numeric_id in group:
            record = result.get(numeric_id)
            if isinstance(record, dict):
                records[f"rs{numeric_id}"] = record
    raw = {
        "requested_rsids": rsids,
        "method": "NCBI dbSNP ESummary",
        "responses": responses,
    }
    return records, raw


def dbsnp_resolution_row(rsid: str, rank: int, record: dict | None) -> dict:
    if not record:
        return {
            "rank": rank, "input_rsid": rsid, "current_rsid": "",
            "was_merged": False, "dbsnp_uid": "", "snp_class": "",
            "spdi": "", "chromosome": "", "chrpos_zero_based": "",
            "merge_text": "", "status": "dbsnp_not_found",
            "error": "No dbSNP ESummary record was returned for this input rsID",
        }
    current_numeric = str(record.get("snp_id", ""))
    if not current_numeric.isdigit():
        current_numeric = rsid.removeprefix("rs")
    current_rsid = f"rs{current_numeric}"
    was_merged = current_rsid != rsid
    return {
        "rank": rank,
        "input_rsid": rsid,
        "current_rsid": current_rsid,
        "was_merged": was_merged,
        "dbsnp_uid": str(record.get("uid", rsid.removeprefix("rs"))),
        "snp_class": record.get("snp_class", ""),
        "spdi": record.get("spdi", ""),
        "chromosome": record.get("chr", ""),
        "chrpos_zero_based": record.get("chrpos", ""),
        "merge_text": record.get("text", ""),
        "status": "merged_to_current_rsid" if was_merged else "current_rsid",
        "error": "",
    }


def fetch_clinvar_links(client: Client, limiter: RateLimiter,
                        current_rsid: str) -> tuple[list[str], dict]:
    """Follow the explicit dbSNP-to-ClinVar ELink relationship."""
    numeric_id = current_rsid.removeprefix("rs")
    response = request_json(client, limiter, "elink.fcgi", {
        "dbfrom": "snp", "db": "clinvar", "id": numeric_id,
        "retmode": "json", "tool": TOOL,
    })
    links = []
    for linkset in response.get("linksets", []):
        if not isinstance(linkset, dict):
            continue
        for linked_database in linkset.get("linksetdbs", []) or []:
            if not isinstance(linked_database, dict):
                continue
            if (str(linked_database.get("dbto", "")).lower() != "clinvar"
                    and linked_database.get("linkname") != "snp_clinvar"):
                continue
            links.extend(str(value) for value in linked_database.get("links", [])
                         if str(value).isdigit())
    return list(dict.fromkeys(links)), response


def fetch_clinvar_summaries(client: Client, limiter: RateLimiter,
                            variation_ids: list[str]) -> tuple[dict[str, dict], dict]:
    """Fetch each unique ClinVar Variation ID once."""
    records: dict[str, dict] = {}
    responses = []
    for group in chunks(variation_ids, SUMMARY_CHUNK_SIZE):
        response = request_json(client, limiter, "esummary.fcgi", {
            "db": "clinvar", "id": ",".join(group), "retmode": "json",
            "version": "2.0", "tool": TOOL,
        })
        responses.append(response)
        result = response.get("result", {})
        for uid in result.get("uids", group):
            record = result.get(str(uid))
            if isinstance(record, dict):
                records[str(uid)] = record
    raw = {
        "requested_variation_ids": variation_ids,
        "method": "NCBI ClinVar ESummary",
        "responses": responses,
    }
    return records, raw


def as_list(value: Any) -> list:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def component_rsids(component: dict) -> list[str]:
    values = []
    for xref in as_list(component.get("variation_xrefs")):
        if not isinstance(xref, dict) or str(xref.get("db_source", "")).lower() != "dbsnp":
            continue
        value = str(xref.get("db_id", "")).removeprefix("rs")
        if value.isdigit() and f"rs{value}" not in values:
            values.append(f"rs{value}")
    return values


def spdi_allele_role(spdi: str) -> str:
    """Distinguish a changed allele from ClinVar's occasional identity allele."""
    try:
        _accession, _position, deleted, inserted = str(spdi).split(":", 3)
    except ValueError:
        return "unavailable"
    return "identity_reference" if deleted.upper() == inserted.upper() else "changed"


def spdi_variant_key(spdi: str) -> str | None:
    """Convert a chromosomal single-base GRCh38 SPDI into the toolkit join key."""
    try:
        accession, position, deleted, inserted = str(spdi).split(":", 3)
        chrom = GRCH38_ACCESSION_TO_CHROM.get(accession.split(".")[0])
        if chrom and len(deleted) == len(inserted) == 1 and deleted != inserted:
            return f"GRCh38:{chrom}:{int(position) + 1}:{deleted.upper()}>{inserted.upper()}"
    except (TypeError, ValueError):
        pass
    return None


def normalize_traits(classification: dict) -> list[dict]:
    traits, seen = [], set()
    for trait in as_list(classification.get("trait_set")):
        if not isinstance(trait, dict):
            continue
        xrefs = [xref for xref in as_list(trait.get("trait_xrefs"))
                 if isinstance(xref, dict)]
        key = (str(trait.get("trait_name", "")), json.dumps(xrefs, sort_keys=True))
        if key in seen:
            continue
        seen.add(key)
        traits.append({"trait_name": key[0], "trait_xrefs": xrefs})
    return traits


def normalize_component(component: dict) -> dict:
    locations = []
    for location in as_list(component.get("variation_loc")):
        if isinstance(location, dict) and location.get("assembly_name") == "GRCh38":
            locations.append({key: location.get(key) for key in
                              ("status", "chr", "start", "stop", "assembly_acc_ver")})
    spdi = str(component.get("canonical_spdi", ""))
    return {
        "measure_id": component.get("measure_id"),
        "variation_name": component.get("variation_name", ""),
        "variant_type": component.get("variant_type", ""),
        "canonical_spdi": spdi,
        "allele_role": spdi_allele_role(spdi),
        "variant_key": spdi_variant_key(spdi),
        "rsids": component_rsids(component),
        "aliases": as_list(component.get("aliases")),
        "grch38_locations": locations,
    }


def normalize_record(rsid: str, current_rsid: str, rank: int, record: dict,
                     clinvar_link_rank: int = 1,
                     allow_linked_set_without_component: bool = False
                     ) -> tuple[dict | None, list[dict]]:
    """Normalize a linked record only when it xrefs the resolved current rsID."""
    components = [item for item in as_list(record.get("variation_set"))
                  if isinstance(item, dict)]
    matching = [component for component in components
                if current_rsid in component_rsids(component)]
    object_type = str(record.get("obj_type", ""))
    explicit_set = any(token in object_type.lower() for token in
                       ("haplotype", "genotype", "compoundheterozygote"))
    if not matching and not (allow_linked_set_without_component and explicit_set):
        return None, []
    is_set = len(components) != 1 or explicit_set
    if not matching:
        scope = "linked_set_without_component_xref"
    else:
        scope = "component_of_set" if is_set else "direct_variant"
    normalized_matching = [normalize_component(component) for component in matching]
    all_component_ids = sorted({item for component in components
                                for item in component_rsids(component)})
    genes = sorted({f"{gene.get('symbol', '')}:{gene.get('geneid', '')}".rstrip(":")
                    for gene in as_list(record.get("genes")) if isinstance(gene, dict)
                    and (gene.get("symbol") or gene.get("geneid"))})
    supporting = record.get("supporting_submissions", {}) or {}
    row: dict[str, Any] = {
        "rank": rank,
        "rsid": rsid,
        "current_rsid": current_rsid,
        "was_merged": current_rsid != rsid,
        "clinvar_link_rank": clinvar_link_rank,
        "clinvar_uid": str(record.get("uid", "")),
        "accession": record.get("accession", ""),
        "accession_version": record.get("accession_version", ""),
        "record_scope": scope,
        "exact_current_dbsnp_xref": bool(matching),
        "matching_allele_roles": (sorted({component["allele_role"]
                                          for component in normalized_matching})
                                  if matching else ["unavailable"]),
        "object_type": object_type,
        "title": record.get("title", ""),
        "variant_keys": sorted({component["variant_key"] for component in normalized_matching
                                if component.get("variant_key")}),
        "matching_components": normalized_matching,
        "component_rsids": all_component_ids,
        "genes": genes,
        "molecular_consequences": as_list(record.get("molecular_consequence_list")),
        "protein_change": record.get("protein_change", ""),
        "n_scv": len(as_list(supporting.get("scv"))),
        "scv_accessions": as_list(supporting.get("scv")),
        "n_rcv": len(as_list(supporting.get("rcv"))),
        "rcv_accessions": as_list(supporting.get("rcv")),
        "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{record.get('uid', '')}/",
    }
    trait_rows = []
    classifications = (
        ("germline", "germline_classification"),
        ("clinical_impact", "clinical_impact_classification"),
        ("oncogenicity", "oncogenicity_classification"),
    )
    for label, key in classifications:
        classification = record.get(key, {}) or {}
        traits = normalize_traits(classification)
        row[f"{label}_classification"] = classification.get("description", "")
        row[f"{label}_review_status"] = classification.get("review_status", "")
        row[f"{label}_last_evaluated"] = classification.get("last_evaluated", "")
        row[f"{label}_traits"] = traits
        for trait in traits:
            trait_rows.append({
                "rsid": rsid,
                "current_rsid": current_rsid,
                "clinvar_uid": row["clinvar_uid"],
                "record_scope": scope,
                "classification_type": label,
                "classification_description": row[f"{label}_classification"],
                "review_status": row[f"{label}_review_status"],
                **trait,
            })
    return row, trait_rows


def aggregate_allele_role(roles: list[str]) -> str:
    unique = sorted(set(roles))
    if not unique:
        return "unavailable"
    return unique[0] if len(unique) == 1 else "mixed"


def sqlite_value(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, bool):
        return int(value)
    return value


def sqlite_columns(fields: list[str], integer_fields: set[str]) -> str:
    return ", ".join(f'"{field}" {"INTEGER" if field in integer_fields else "TEXT"}'
                     for field in fields)


def insert_rows(database: sqlite3.Connection, table: str, fields: list[str],
                rows: list[dict]) -> None:
    quoted = ", ".join(f'"{field}"' for field in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{table}" ({quoted}) VALUES ({placeholders})',
        [tuple(sqlite_value(row.get(field, "")) for field in fields) for row in rows],
    )


def write_sqlite(path: Path, dbsnp_rows: list[dict], records: list[dict],
                 links: list[dict], statuses: list[dict], traits: list[dict],
                 manifest: dict) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    try:
        with database:
            database.execute(
                f"CREATE TABLE dbsnp_resolution ({sqlite_columns(DBSNP_FIELDS, {'rank', 'was_merged'})}, "
                "PRIMARY KEY (input_rsid))"
            )
            database.execute(
                f"CREATE TABLE records ({sqlite_columns(RECORD_FIELDS, {'rank', 'was_merged', 'clinvar_link_rank', 'exact_current_dbsnp_xref', 'n_scv', 'n_rcv'})}, "
                "PRIMARY KEY (rsid, clinvar_uid))"
            )
            database.execute(
                f"CREATE TABLE links ({sqlite_columns(LINK_FIELDS, {'rank', 'was_merged', 'clinvar_link_rank', 'matching_component_count', 'exact_current_dbsnp_xref'})}, "
                "PRIMARY KEY (rsid, clinvar_uid))"
            )
            status_integer = {"rank", "was_merged"} | {field for field in STATUS_FIELDS
                                                        if field.startswith("n_")}
            database.execute(
                f"CREATE TABLE panel_status ({sqlite_columns(STATUS_FIELDS, status_integer)}, "
                "PRIMARY KEY (rsid))"
            )
            database.execute(f"CREATE TABLE traits ({sqlite_columns(TRAIT_FIELDS, set())})")
            database.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)")
            insert_rows(database, "dbsnp_resolution", DBSNP_FIELDS, dbsnp_rows)
            insert_rows(database, "records", RECORD_FIELDS, records)
            insert_rows(database, "links", LINK_FIELDS, links)
            insert_rows(database, "panel_status", STATUS_FIELDS, statuses)
            insert_rows(database, "traits", TRAIT_FIELDS, traits)
            database.execute("INSERT INTO metadata VALUES (?, ?)",
                             ("manifest", json.dumps(manifest, ensure_ascii=False,
                                                      sort_keys=True)))
            integrity = database.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise ValueError(f"SQLite integrity check failed: {integrity}")
    finally:
        database.close()
    partial.replace(path)


def write_readme(out: Path, manifest: dict) -> None:
    text = f"""# ClinVar dataset for the QIAGEN 20-rsID panel

This folder contains an NCBI dbSNP-to-ClinVar snapshot for the exact 20 public
test rsIDs. It contains no patient genotypes.

- Input rsIDs: {manifest['n_input_rsids']}
- Merged input rsIDs resolved by dbSNP: {manifest['n_merged_input_rsids']}
- Explicit rsID-to-ClinVar links: {manifest['n_clinvar_links']}
- Direct ClinVar records: {manifest['n_direct_records']}
- Related haplotype/genotype records: {manifest['n_related_set_records']}
- Linked sets without returned component xrefs: {manifest['n_linked_sets_without_component_xref']}
- Input rsIDs with at least one direct record: {manifest['n_rsids_with_direct_records']}

## Files

`dbsnp_resolution.tsv` maps each supplied ID to its current dbSNP ID. The
supplied ID remains the primary panel key. `rsid_clinvar_links.tsv` is the
auditable link table from each panel ID to every linked ClinVar Variation ID.
`clinvar_records.tsv` contains direct single-variant records.
`related_set_records.tsv` contains multi-variant records in which the rsID is
only one component; do not assign their set-level classifications to that SNP
alone. `panel_status.tsv` keeps all 20 inputs in original order. JSON copies
preserve nested values, `raw/` preserves API responses, and `clinvar.sqlite`
contains `dbsnp_resolution`, `links`, `records`, `panel_status`, `traits`, and
`metadata` tables.

`traits.tsv` expands the trait labels returned in each VCV-level ESummary
record. Its classification fields are VCV aggregates repeated per trait; they
are not condition-specific RCV assertion records.

`allele_role=identity_reference` marks a linked ClinVar record whose canonical
SPDI has identical deleted and inserted sequence. Such a record is retained
for provenance but is not a changed alternate allele. `variant_keys` are only
created for changed single-base GRCh38 alleles.

Some explicit ELink relationships lead to a haplotype/genotype record for
which ClinVar ESummary does not return component objects. Those records remain
in `related_set_records` with `record_scope=linked_set_without_component_xref`
and `exact_current_dbsnp_xref=false`; this preserves the link without treating
its set-level classification as an allele-level interpretation.

ClinVar reports submitted interpretations and review status. Germline,
clinical-impact, and oncogenicity classifications remain separate. These data
are not a diagnosis, a personal risk estimate, or an independently verified
truth set. Interpret every classification with its exact allele, condition,
review status, and evaluation date.

Sources: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/ and
https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def empty_status(rank: int, rsid: str, resolution: dict, status: str,
                 error: str = "") -> dict:
    return {
        "rank": rank,
        "rsid": rsid,
        "current_rsid": resolution.get("current_rsid", ""),
        "was_merged": resolution.get("was_merged", False),
        "dbsnp_status": resolution.get("status", ""),
        "dbsnp_snp_class": resolution.get("snp_class", ""),
        "dbsnp_spdi": resolution.get("spdi", ""),
        "status": status,
        "n_clinvar_links": 0,
        "n_linked_summaries": 0,
        "n_exact_records": 0,
        "n_direct_records": 0,
        "n_related_set_records": 0,
        "n_rejected_linked_records": 0,
        "direct_accessions": [],
        "direct_germline_classifications": [],
        "direct_review_statuses": [],
        "error": error,
    }


def run(args: argparse.Namespace, client: Client | None = None) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    if args.max_records_per_rsid < 1:
        raise ValueError("--max-records-per-rsid must be positive")
    rsids = select_panel_rsids(args.input, require_exact=True)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    client = client or Client(args.cache, args.offline, args.refresh)
    limiter = RateLimiter(0 if args.offline else 0.34)

    print(f"Resolving {len(rsids)} supplied IDs with NCBI dbSNP...", flush=True)
    dbsnp_records, dbsnp_raw = fetch_dbsnp_summaries(client, limiter, rsids)
    write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)
    dbsnp_rows = [dbsnp_resolution_row(rsid, rank, dbsnp_records.get(rsid))
                   for rank, rsid in enumerate(rsids, 1)]
    resolution_by_input = {row["input_rsid"]: row for row in dbsnp_rows}

    current_rsids = list(dict.fromkeys(row["current_rsid"] for row in dbsnp_rows
                                       if row["current_rsid"]))
    links_by_current: dict[str, list[str]] = {}
    link_errors: dict[str, str] = {}
    for index, current_rsid in enumerate(current_rsids, 1):
        print(f"[{index}/{len(current_rsids)}] Following dbSNP ClinVar links for "
              f"{current_rsid}...", flush=True)
        try:
            links, raw = fetch_clinvar_links(client, limiter, current_rsid)
            if len(links) > args.max_records_per_rsid:
                raise ValueError(
                    f"{current_rsid} has {len(links)} ClinVar links, exceeding "
                    f"--max-records-per-rsid={args.max_records_per_rsid}"
                )
            links_by_current[current_rsid] = links
            write_json(out / "raw" / "elink" / f"{current_rsid}.json", raw)
        except Exception as exc:
            link_errors[current_rsid] = str(exc)
            links_by_current[current_rsid] = []

    unique_variation_ids = list(dict.fromkeys(
        uid for current_rsid in current_rsids for uid in links_by_current[current_rsid]
    ))
    print(f"Fetching {len(unique_variation_ids)} unique ClinVar summaries...", flush=True)
    clinvar_records, clinvar_raw = fetch_clinvar_summaries(
        client, limiter, unique_variation_ids
    )
    write_json(out / "raw" / "clinvar_esummary.json", clinvar_raw)

    direct_records: list[dict] = []
    related_records: list[dict] = []
    link_rows: list[dict] = []
    statuses: list[dict] = []
    trait_rows: list[dict] = []
    errors: list[dict] = []

    for rank, rsid in enumerate(rsids, 1):
        resolution = resolution_by_input[rsid]
        current_rsid = resolution["current_rsid"]
        if not current_rsid:
            error = resolution["error"]
            errors.append({"rank": rank, "rsid": rsid,
                           "stage": "dbsnp_resolution", "error": error})
            statuses.append(empty_status(rank, rsid, resolution, "dbsnp_not_found", error))
            continue
        if current_rsid in link_errors:
            error = link_errors[current_rsid]
            errors.append({"rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                           "stage": "dbsnp_clinvar_link", "error": error})
            statuses.append(empty_status(rank, rsid, resolution, "lookup_failed", error))
            continue

        variation_ids = links_by_current[current_rsid]
        local_direct, local_related, local_traits = [], [], []
        rejected = 0
        linked_summaries = 0
        for link_rank, uid in enumerate(variation_ids, 1):
            record = clinvar_records.get(uid)
            if record is None:
                rejected += 1
                link_rows.append({
                    "rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                    "was_merged": current_rsid != rsid,
                    "clinvar_link_rank": link_rank, "clinvar_uid": uid,
                    "accession": "", "accession_version": "",
                    "record_scope": "summary_missing", "allele_role": "unavailable",
                    "matching_variant_keys": [], "matching_component_count": 0,
                    "exact_current_dbsnp_xref": False,
                    "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/",
                })
                continue
            linked_summaries += 1
            row, rows_for_traits = normalize_record(
                rsid, current_rsid, rank, record, link_rank,
                allow_linked_set_without_component=True,
            )
            if row is None:
                rejected += 1
                link_rows.append({
                    "rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                    "was_merged": current_rsid != rsid,
                    "clinvar_link_rank": link_rank, "clinvar_uid": uid,
                    "accession": record.get("accession", ""),
                    "accession_version": record.get("accession_version", ""),
                    "record_scope": "linked_without_exact_current_dbsnp_xref",
                    "allele_role": "unavailable", "matching_variant_keys": [],
                    "matching_component_count": 0,
                    "exact_current_dbsnp_xref": False,
                    "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/",
                })
                continue
            local_traits.extend(rows_for_traits)
            if row["record_scope"] == "direct_variant":
                local_direct.append(row)
            else:
                local_related.append(row)
            link_rows.append({
                "rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                "was_merged": current_rsid != rsid,
                "clinvar_link_rank": link_rank, "clinvar_uid": uid,
                "accession": row["accession"],
                "accession_version": row["accession_version"],
                "record_scope": row["record_scope"],
                "allele_role": aggregate_allele_role(row["matching_allele_roles"]),
                "matching_variant_keys": row["variant_keys"],
                "matching_component_count": len(row["matching_components"]),
                "exact_current_dbsnp_xref": row["exact_current_dbsnp_xref"],
                "source_url": row["source_url"],
            })

        direct_records.extend(local_direct)
        related_records.extend(local_related)
        trait_rows.extend(local_traits)
        if local_direct:
            status = "direct_records_found"
        elif local_related:
            status = "related_set_records_only"
        elif variation_ids:
            status = "linked_records_without_exact_dbsnp_xref"
        else:
            status = "no_clinvar_record"
        status_row = empty_status(rank, rsid, resolution, status)
        status_row.update({
            "n_clinvar_links": len(variation_ids),
            "n_linked_summaries": linked_summaries,
            "n_exact_records": sum(row["exact_current_dbsnp_xref"]
                                   for row in local_direct + local_related),
            "n_direct_records": len(local_direct),
            "n_related_set_records": len(local_related),
            "n_rejected_linked_records": rejected,
            "direct_accessions": [row["accession_version"] or row["accession"]
                                  for row in local_direct],
            "direct_germline_classifications": sorted({row["germline_classification"]
                                                        for row in local_direct
                                                        if row["germline_classification"]}),
            "direct_review_statuses": sorted({row["germline_review_status"]
                                               for row in local_direct
                                               if row["germline_review_status"]}),
        })
        statuses.append(status_row)

    all_records = direct_records + related_records
    input_path = Path(args.input)
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "source_databases": ["NCBI dbSNP", "NCBI ClinVar"],
        "api": "NCBI E-utilities dbSNP ESummary -> ELink -> ClinVar ESummary",
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "n_input_rsids": len(rsids),
        "n_merged_input_rsids": sum(bool(row["was_merged"]) for row in dbsnp_rows),
        "n_clinvar_links": len(link_rows),
        "n_unique_clinvar_records": len(clinvar_records),
        "n_rsids_with_direct_records": sum(row["n_direct_records"] > 0 for row in statuses),
        "n_direct_records": len(direct_records),
        "n_related_set_records": len(related_records),
        "n_linked_sets_without_component_xref": sum(
            row["record_scope"] == "linked_set_without_component_xref"
            for row in related_records
        ),
        "n_trait_rows": len(trait_rows),
        "n_lookup_errors": len(errors),
        "source_provenance": list(client.manifest.values()),
        "notes": [
            "Supplied rsIDs are preserved while dbSNP merged/current IDs are recorded separately",
            "Only explicit dbSNP-to-ClinVar links are downloaded",
            "Exact component xrefs are verified; explicitly linked sets without returned components are labeled unverified",
            "Direct records and multi-variant set records are deliberately separated",
            "Identity/reference alleles are retained and explicitly marked",
            "Germline, clinical-impact, and oncogenicity classifications remain separate",
            "Trait rows repeat VCV-level aggregate classifications and are not RCV-level assertions",
            "The input contains public loci, not an individual's genotype",
        ],
    }

    write_json(out / "dbsnp_resolution.json", dbsnp_rows)
    write_tsv(out / "dbsnp_resolution.tsv", dbsnp_rows, DBSNP_FIELDS)
    write_json(out / "rsid_clinvar_links.json", link_rows)
    write_tsv(out / "rsid_clinvar_links.tsv", link_rows, LINK_FIELDS)
    write_json(out / "clinvar_records.json", direct_records)
    write_tsv(out / "clinvar_records.tsv", direct_records, RECORD_FIELDS)
    write_json(out / "related_set_records.json", related_records)
    write_tsv(out / "related_set_records.tsv", related_records, RECORD_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "traits.json", trait_rows)
    write_tsv(out / "traits.tsv", trait_rows, TRAIT_FIELDS)
    write_json(out / "errors.json", errors)
    write_json(out / "run_manifest.json", manifest)
    write_sqlite(out / "clinvar.sqlite", dbsnp_rows, all_records, link_rows,
                 statuses, trait_rows, manifest)
    write_readme(out, manifest)
    print(f"ClinVar dataset complete: {len(rsids)} input rsIDs, {len(link_rows)} links, "
          f"{len(direct_records)} direct records, {len(related_records)} related set "
          f"records. See {out}")
    return 1 if errors else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="Exact 20-rsID QIAGEN TSV in the original order")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument("--cache", default=".cache", help="Shared public-API response cache")
    parser.add_argument("--offline", action="store_true", help="Use cache only")
    parser.add_argument("--refresh", action="store_true", help="Refresh cached requests")
    parser.add_argument("--max-records-per-rsid", type=int, default=1000,
                        help="Safety limit; fail an rsID rather than silently truncate")
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
