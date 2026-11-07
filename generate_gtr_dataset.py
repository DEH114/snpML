"""Build an NCBI Genetic Testing Registry dataset for the 20 QIAGEN rsIDs.

GTR is a registry of laboratory tests, not a variant-classification database.
This generator resolves current dbSNP identifiers, performs exact quoted-word
Entrez searches for each supplied/current rsID, and downloads deduplicated GTR
test summaries. A match means the registered test is indexed for that rsID; it
does not establish clinical validity, FDA review, or a patient result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
from pathlib import Path
from typing import Any

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import Client, utc_now, write_json, write_tsv
from generate_clinvar_dataset import (
    DBSNP_FIELDS,
    RateLimiter,
    chunks,
    dbsnp_resolution_row,
    fetch_dbsnp_summaries,
    request_json,
    sqlite_value,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "gtr_dataset"
SEARCH_PAGE_SIZE = 100
SUMMARY_CHUNK_SIZE = 100

TEST_FIELDS = [
    "gtr_uid", "accession", "test_name", "test_type", "source", "offerer",
    "offerer_id", "offerer_location", "conditions", "analytes", "methods",
    "analytical_validity", "clinical_validity", "target_population",
    "certifications", "test_targets", "variants", "test_purposes",
    "test_development", "specimens", "country", "test_url", "order_url",
    "specimen_url", "test_code", "cpt_code", "summary", "source_url",
]
LINK_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "gtr_uid", "accession",
    "test_name", "test_type", "offerer", "matched_query_rsids",
    "relation_basis", "source_url",
]
STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "status",
    "n_query_terms", "n_search_hits", "n_test_summaries", "accessions",
    "matched_query_rsids", "error",
]
CONDITION_FIELDS = [
    "gtr_uid", "accession", "condition_name", "condition_acronym",
    "condition_cui",
]
VARIANT_FIELDS = ["gtr_uid", "accession", "target_type", "target_value"]


def as_list(value: Any) -> list:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def fetch_gtr_ids(client: Client, limiter: RateLimiter, query_rsid: str,
                  max_records: int = 1000) -> tuple[list[str], dict]:
    """Search the exact quoted rsID token in GTR's Entrez WORD index."""
    ids: list[str] = []
    pages = []
    expected_count: int | None = None
    retstart = 0
    term = f'"{query_rsid}"[WORD]'
    while expected_count is None or retstart < expected_count:
        page = request_json(client, limiter, "esearch.fcgi", {
            "db": "gtr", "term": term, "retmode": "json",
            "retstart": retstart, "retmax": SEARCH_PAGE_SIZE,
            "tool": "BioHack_SNP_Toolkit",
        })
        pages.append(page)
        result = page.get("esearchresult", {})
        translation = str(result.get("querytranslation", "")).lower()
        if query_rsid.lower() not in translation or "word" not in translation:
            raise ValueError(
                f"GTR changed or broadened the exact query translation for {query_rsid}: "
                f"{result.get('querytranslation', '')!r}"
            )
        try:
            count = int(result.get("count", 0))
        except (TypeError, ValueError) as exc:
            raise ValueError("GTR ESearch returned an invalid count") from exc
        if expected_count is None:
            expected_count = count
            if count > max_records:
                raise ValueError(
                    f"{query_rsid} returned {count} GTR tests, exceeding "
                    f"--max-records-per-rsid={max_records}"
                )
        page_ids = [str(value) for value in result.get("idlist", [])
                    if str(value).isdigit()]
        ids.extend(page_ids)
        if not page_ids:
            break
        retstart += len(page_ids)
    unique_ids = list(dict.fromkeys(ids))
    return unique_ids, {
        "query_rsid": query_rsid,
        "query": term,
        "expected_count": expected_count or 0,
        "returned_uids": unique_ids,
        "pages": pages,
    }


def fetch_gtr_summaries(client: Client, limiter: RateLimiter,
                        uids: list[str]) -> tuple[dict[str, dict], dict]:
    records: dict[str, dict] = {}
    responses = []
    for group in chunks(uids, SUMMARY_CHUNK_SIZE):
        response = request_json(client, limiter, "esummary.fcgi", {
            "db": "gtr", "id": ",".join(group), "retmode": "json",
            "version": "2.0", "tool": "BioHack_SNP_Toolkit",
        })
        responses.append(response)
        result = response.get("result", {})
        for uid in result.get("uids", group):
            record = result.get(str(uid))
            if isinstance(record, dict):
                records[str(uid)] = record
    return records, {
        "requested_uids": uids,
        "method": "NCBI GTR ESummary",
        "responses": responses,
    }


def normalize_test(record: dict) -> dict:
    uid = str(record.get("uid", record.get("id", "")))
    accession = str(record.get("accession", ""))
    return {
        "gtr_uid": uid,
        "accession": accession,
        "test_name": record.get("testname", ""),
        "test_type": record.get("testtype", ""),
        "source": record.get("source", ""),
        "offerer": record.get("offerer", ""),
        "offerer_id": record.get("offererid", ""),
        "offerer_location": record.get("offererlocation", {}) or {},
        "conditions": as_list(record.get("conditionlist") or record.get("conditionlist2")),
        "analytes": as_list(record.get("analytes")) or as_list(record.get("genelist")),
        "methods": as_list(record.get("method")),
        "analytical_validity": record.get("analyticalvalidity", ""),
        "clinical_validity": record.get("clinicalvalidity", {}) or {},
        "target_population": record.get("targetpopulation", "")
                             or record.get("targetpopulation2", ""),
        "certifications": as_list(record.get("certifications")),
        "test_targets": as_list(record.get("testtargetlist"))
                        or as_list(record.get("testtargetlist2")),
        "variants": as_list(record.get("variant")),
        "test_purposes": as_list(record.get("testpurpose")),
        "test_development": record.get("testdevelopment", ""),
        "specimens": as_list(record.get("specimens")),
        "country": record.get("country", ""),
        "test_url": record.get("testurl", ""),
        "order_url": record.get("orderurl", ""),
        "specimen_url": record.get("specimenurl", ""),
        "test_code": record.get("testcode", ""),
        "cpt_code": record.get("cptcode", ""),
        "summary": record.get("summary", ""),
        "source_url": f"https://www.ncbi.nlm.nih.gov/gtr/tests/{accession}/"
                      if accession else f"https://www.ncbi.nlm.nih.gov/gtr/?term={uid}",
    }


def condition_rows(tests: list[dict]) -> list[dict]:
    rows = []
    for test in tests:
        seen = set()
        for condition in test["conditions"]:
            if not isinstance(condition, dict):
                condition = {"name": str(condition)}
            key = (str(condition.get("name", "")), str(condition.get("cui", "")))
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "gtr_uid": test["gtr_uid"], "accession": test["accession"],
                "condition_name": condition.get("name", ""),
                "condition_acronym": condition.get("acronym", ""),
                "condition_cui": condition.get("cui", ""),
            })
    return rows


def variant_target_rows(tests: list[dict]) -> list[dict]:
    rows = []
    for test in tests:
        seen = set()
        for target_type, values in (("registered_target", test["test_targets"]),
                                    ("variant_hgvs", test["variants"])):
            for value in values:
                text = value if isinstance(value, str) else json.dumps(value, sort_keys=True)
                key = (target_type, text)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    "gtr_uid": test["gtr_uid"], "accession": test["accession"],
                    "target_type": target_type, "target_value": text,
                })
    return rows


def sqlite_columns(fields: list[str], integers: set[str]) -> str:
    return ", ".join(f'"{field}" {"INTEGER" if field in integers else "TEXT"}'
                     for field in fields)


def insert_rows(database: sqlite3.Connection, table: str, fields: list[str],
                rows: list[dict]) -> None:
    quoted = ", ".join(f'"{field}"' for field in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{table}" ({quoted}) VALUES ({placeholders})',
        [tuple(sqlite_value(row.get(field, "")) for field in fields) for row in rows],
    )


def write_sqlite(path: Path, dbsnp_rows: list[dict], tests: list[dict],
                 links: list[dict], statuses: list[dict], conditions: list[dict],
                 targets: list[dict], manifest: dict) -> None:
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
                f"CREATE TABLE tests ({sqlite_columns(TEST_FIELDS, {'offerer_id'})}, "
                "PRIMARY KEY (gtr_uid))"
            )
            database.execute(
                f"CREATE TABLE links ({sqlite_columns(LINK_FIELDS, {'rank', 'was_merged'})}, "
                "PRIMARY KEY (rsid, gtr_uid))"
            )
            status_ints = {"rank", "was_merged", "n_query_terms", "n_search_hits",
                           "n_test_summaries"}
            database.execute(
                f"CREATE TABLE panel_status ({sqlite_columns(STATUS_FIELDS, status_ints)}, "
                "PRIMARY KEY (rsid))"
            )
            database.execute(f"CREATE TABLE conditions ({sqlite_columns(CONDITION_FIELDS, set())})")
            database.execute(f"CREATE TABLE variant_targets ({sqlite_columns(VARIANT_FIELDS, set())})")
            database.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)")
            insert_rows(database, "dbsnp_resolution", DBSNP_FIELDS, dbsnp_rows)
            insert_rows(database, "tests", TEST_FIELDS, tests)
            insert_rows(database, "links", LINK_FIELDS, links)
            insert_rows(database, "panel_status", STATUS_FIELDS, statuses)
            insert_rows(database, "conditions", CONDITION_FIELDS, conditions)
            insert_rows(database, "variant_targets", VARIANT_FIELDS, targets)
            database.execute("INSERT INTO metadata VALUES (?, ?)",
                             ("manifest", json.dumps(manifest, sort_keys=True)))
            integrity = database.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise ValueError(f"SQLite integrity check failed: {integrity}")
    finally:
        database.close()
    partial.replace(path)


def write_readme(out: Path, manifest: dict) -> None:
    text = f"""# GTR registered-test matches for the QIAGEN 20-rsID panel

This snapshot searches NCBI's Genetic Testing Registry (GTR) for the exact 20
public rsIDs. It contains registered laboratory test summaries, not patient
records or variant classifications.

- Input rsIDs: {manifest['n_input_rsids']}
- rsIDs with at least one indexed test match: {manifest['n_rsids_with_matches']}
- rsID-to-test relationships: {manifest['n_test_links']}
- Unique registered tests: {manifest['n_unique_tests']}

`panel_status.tsv` contains all 20 inputs in order. `gtr_test_matches.tsv`
preserves which exact quoted rsID search returned each test. `gtr_tests.tsv`
contains deduplicated test summaries; `conditions.tsv` and
`variant_targets.tsv` expand selected nested fields. JSON copies preserve
nested values, `raw/` preserves Entrez responses, and `gtr.sqlite` contains the
normalized tables.

GTR does not expose a dedicated rsID field or dbSNP-to-GTR link. Therefore the
relationship basis is an exact quoted `[WORD]` Entrez index match. An rsID may
not appear literally in ESummary even when the search index links the test.
The test may assay a gene, panel, or equivalent HGVS target rather than only
that SNP.

GTR content is submitted by laboratories and is not independently verified by
NIH. A match is not a pathogenicity assertion, FDA approval, recommendation,
or proof of clinical validity. Absence of a match is not evidence that testing
is unavailable.

Official access documentation:
https://www.ncbi.nlm.nih.gov/gtr/docs/maintenance_use/
"""
    (out / "README.md").write_text(text, encoding="utf-8")


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
    dbsnp_rows = [dbsnp_resolution_row(rsid, rank, dbsnp_records.get(rsid))
                  for rank, rsid in enumerate(rsids, 1)]
    resolution_by_input = {row["input_rsid"]: row for row in dbsnp_rows}
    write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)

    matches_by_input: dict[str, dict[str, list[str]]] = {}
    search_errors: dict[str, list[str]] = {}
    searched: dict[str, tuple[list[str], dict] | Exception] = {}
    query_plan = []
    for rsid in rsids:
        current = resolution_by_input[rsid]["current_rsid"]
        query_ids = list(dict.fromkeys([value for value in (rsid, current) if value]))
        query_plan.append((rsid, query_ids))
        for query_rsid in query_ids:
            if query_rsid in searched:
                continue
            print(f"Searching GTR for exact token {query_rsid}...", flush=True)
            try:
                ids, raw = fetch_gtr_ids(client, limiter, query_rsid,
                                         args.max_records_per_rsid)
                searched[query_rsid] = (ids, raw)
                write_json(out / "raw" / "esearch" / f"{query_rsid}.json", raw)
            except Exception as exc:
                searched[query_rsid] = exc

    for rsid, query_ids in query_plan:
        by_uid: dict[str, list[str]] = {}
        local_errors = []
        for query_rsid in query_ids:
            result = searched[query_rsid]
            if isinstance(result, Exception):
                local_errors.append(f"{query_rsid}: {result}")
                continue
            ids, _raw = result
            for uid in ids:
                by_uid.setdefault(uid, []).append(query_rsid)
        matches_by_input[rsid] = by_uid
        search_errors[rsid] = local_errors

    unique_uids = list(dict.fromkeys(uid for rsid in rsids
                                     for uid in matches_by_input[rsid]))
    print(f"Fetching {len(unique_uids)} unique GTR test summaries...", flush=True)
    summary_records, summary_raw = fetch_gtr_summaries(client, limiter, unique_uids)
    write_json(out / "raw" / "esummary.json", summary_raw)
    tests = [normalize_test(summary_records[uid]) for uid in unique_uids
             if uid in summary_records]
    tests_by_uid = {row["gtr_uid"]: row for row in tests}

    links, statuses, errors = [], [], []
    for rank, rsid in enumerate(rsids, 1):
        resolution = resolution_by_input[rsid]
        current = resolution["current_rsid"]
        by_uid = matches_by_input[rsid]
        local_links = []
        missing_summaries = []
        for uid, query_rsids in by_uid.items():
            test = tests_by_uid.get(uid)
            if not test:
                missing_summaries.append(uid)
                continue
            if rsid in query_rsids and current != rsid and current in query_rsids:
                basis = "entrez_exact_word_input_and_current_rsid"
            elif current in query_rsids and current != rsid:
                basis = "entrez_exact_word_current_rsid"
            else:
                basis = "entrez_exact_word_input_rsid"
            local_links.append({
                "rank": rank, "rsid": rsid, "current_rsid": current,
                "was_merged": resolution["was_merged"], "gtr_uid": uid,
                "accession": test["accession"], "test_name": test["test_name"],
                "test_type": test["test_type"], "offerer": test["offerer"],
                "matched_query_rsids": query_rsids, "relation_basis": basis,
                "source_url": test["source_url"],
            })
        links.extend(local_links)
        local_errors = list(search_errors[rsid])
        if missing_summaries:
            local_errors.append("Missing GTR summaries for UIDs: " + ",".join(missing_summaries))
        if local_errors:
            errors.append({"rank": rank, "rsid": rsid, "stage": "gtr_lookup",
                           "error": "; ".join(local_errors)})
        if local_links and local_errors:
            status = "matches_found_with_partial_error"
        elif local_links:
            status = "registered_tests_found"
        elif local_errors:
            status = "lookup_failed"
        else:
            status = "no_gtr_test_match"
        query_ids = next(values for panel_id, values in query_plan if panel_id == rsid)
        statuses.append({
            "rank": rank, "rsid": rsid, "current_rsid": current,
            "was_merged": resolution["was_merged"], "status": status,
            "n_query_terms": len(query_ids), "n_search_hits": len(by_uid),
            "n_test_summaries": len(local_links),
            "accessions": [row["accession"] for row in local_links],
            "matched_query_rsids": sorted({query for row in local_links
                                            for query in row["matched_query_rsids"]}),
            "error": "; ".join(local_errors),
        })

    conditions = condition_rows(tests)
    targets = variant_target_rows(tests)
    input_path = Path(args.input)
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "source_databases": ["NCBI dbSNP", "NCBI Genetic Testing Registry"],
        "api": "NCBI E-utilities dbSNP ESummary + GTR exact WORD ESearch/ESummary",
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "n_input_rsids": len(rsids),
        "n_rsids_with_matches": sum(row["n_test_summaries"] > 0 for row in statuses),
        "n_test_links": len(links),
        "n_unique_tests": len(tests),
        "n_conditions": len(conditions),
        "n_variant_targets": len(targets),
        "n_lookup_errors": len(errors),
        "source_provenance": list(client.manifest.values()),
        "notes": [
            "GTR matches are registered laboratory tests, not variant evidence or patient results",
            "Relationship provenance is an exact quoted Entrez WORD-index match",
            "Laboratory submissions are not independently verified by NIH",
            "A missing GTR match is not evidence that testing is unavailable",
        ],
    }

    write_json(out / "dbsnp_resolution.json", dbsnp_rows)
    write_tsv(out / "dbsnp_resolution.tsv", dbsnp_rows, DBSNP_FIELDS)
    write_json(out / "gtr_tests.json", tests)
    write_tsv(out / "gtr_tests.tsv", tests, TEST_FIELDS)
    write_json(out / "gtr_test_matches.json", links)
    write_tsv(out / "gtr_test_matches.tsv", links, LINK_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "conditions.json", conditions)
    write_tsv(out / "conditions.tsv", conditions, CONDITION_FIELDS)
    write_json(out / "variant_targets.json", targets)
    write_tsv(out / "variant_targets.tsv", targets, VARIANT_FIELDS)
    write_json(out / "errors.json", errors)
    write_json(out / "run_manifest.json", manifest)
    write_sqlite(out / "gtr.sqlite", dbsnp_rows, tests, links, statuses,
                 conditions, targets, manifest)
    write_readme(out, manifest)
    print(f"GTR dataset complete: {len(rsids)} rsIDs, {len(links)} test links, "
          f"{len(tests)} unique tests. See {out}")
    return 1 if errors else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="Exact 20-rsID QIAGEN TSV in original order")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument("--cache", default=".cache", help="Shared public-API response cache")
    parser.add_argument("--offline", action="store_true", help="Use cache only")
    parser.add_argument("--refresh", action="store_true", help="Refresh cached requests")
    parser.add_argument("--max-records-per-rsid", type=int, default=1000)
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
