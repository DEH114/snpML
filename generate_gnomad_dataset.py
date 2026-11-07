"""Build a gnomAD frequency dataset for the exact 20 QIAGEN challenge rsIDs.

The supplied IDs are first resolved to current dbSNP identifiers. The current
IDs are searched in the official gnomAD v4 GraphQL API, because one rsID may
represent several alleles and indels may be normalized differently by gnomAD.
Each returned gnomAD variant ID is then fetched explicitly and its dbSNP xref
is verified before frequencies are exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import requests

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import APIError, Client, utc_now, write_json, write_tsv
from generate_clinvar_dataset import (DBSNP_FIELDS, RateLimiter,
                                      dbsnp_resolution_row,
                                      fetch_dbsnp_summaries)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / "results" / "qiagen_20_snp_dataset" / "gnomad_dataset"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
DETAIL_CHUNK_SIZE = 8
RETRYABLE_HTTP = {429, 500, 502, 503, 504}
VARIANT_ID_RE = re.compile(
    r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y)-[1-9][0-9]*-[ACGT]+-[ACGT]+$"
)

LINK_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "gnomad_variant_rank",
    "variant_id", "detail_available", "exact_current_rsid_xref", "status",
    "source_url",
]
VARIANT_FIELDS = [
    "variant_id", "reference_genome", "chrom", "pos", "ref", "alt", "caid",
    "rsids", "flags", "source_url",
]
CALLSET_FIELDS = [
    "variant_id", "callset", "ac", "an", "af", "reported_af",
    "homozygote_count", "hemizygote_count", "filters", "flags",
    "faf95_popmax_population", "faf95_popmax",
]
POPULATION_FIELDS = [
    "variant_id", "callset", "population_id", "population_scope", "ac", "an",
    "af", "homozygote_count", "hemizygote_count",
]
STATUS_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged", "dbsnp_status",
    "dbsnp_snp_class", "dbsnp_spdi", "status", "n_search_hits",
    "n_verified_variants", "n_detail_missing", "n_verification_failed",
    "gnomad_variant_ids", "error",
]


class GnomadGraphQLClient:
    """Cached POST-only client for the public gnomAD GraphQL endpoint.

    Cache keys include the complete canonical JSON body. HTTP throttling and
    transient GraphQL rate-limit errors are retried; ordinary GraphQL errors
    are returned so callers can retain partial data and an audit trail.
    """

    def __init__(self, cache: str | Path, offline: bool = False,
                 refresh: bool = False, timeout: int = 60,
                 session: requests.Session | None = None,
                 sleep: Any = time.sleep):
        if offline and refresh:
            raise ValueError("--offline and --refresh are incompatible")
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.refresh = refresh
        self.timeout = timeout
        self.session = session or requests.Session()
        self.session.headers.update({
            "User-Agent": "BioHack-SNP-Toolkit/1.0 (research prototype)",
            "Accept": "application/json",
        })
        self.sleep = sleep
        self.manifest: dict[str, dict] = {}

    @staticmethod
    def _is_transient_graphql_error(data: dict) -> bool:
        return any(
            "rate limit exceeded" in str(error.get("message", "")).lower()
            for error in data.get("errors", []) if isinstance(error, dict)
        )

    def post(self, query: str, variables: dict[str, Any] | None = None,
             operation_name: str | None = None) -> dict:
        body: dict[str, Any] = {"query": query, "variables": variables or {}}
        if operation_name:
            body["operationName"] = operation_name
        canonical = json.dumps(body, sort_keys=True, separators=(",", ":"),
                               ensure_ascii=False)
        key = hashlib.sha256((GNOMAD_API + "\n" + canonical).encode()).hexdigest()
        payload_path = self.cache / f"{key}.response.json"
        metadata_path = self.cache / f"{key}.meta.json"
        if payload_path.exists() and metadata_path.exists() and not self.refresh:
            data = json.loads(payload_path.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise APIError(f"Invalid cached gnomAD response: {payload_path}")
            self.manifest[key] = json.loads(metadata_path.read_text(encoding="utf-8"))
            return data
        if self.offline:
            raise APIError(
                f"Offline gnomAD cache miss for operation {operation_name or 'anonymous'}"
            )

        for attempt in range(4):
            try:
                response = self.session.post(
                    GNOMAD_API, json=body,
                    headers={"Content-Type": "application/json"},
                    timeout=(10, self.timeout),
                )
            except requests.RequestException as exc:
                if attempt == 3:
                    raise APIError(f"gnomAD network error: {exc}") from exc
                self.sleep(2 ** attempt)
                continue

            if response.status_code in RETRYABLE_HTTP and attempt < 3:
                try:
                    wait = min(float(response.headers.get("Retry-After", 2 ** attempt)), 60)
                except (TypeError, ValueError):
                    wait = 2 ** attempt
                self.sleep(wait)
                continue
            if not response.ok:
                raise APIError(
                    f"gnomAD HTTP {response.status_code}: {response.text[:350]}",
                    response.status_code,
                )
            try:
                data = json.loads(response.text)
            except json.JSONDecodeError as exc:
                if attempt < 3:
                    self.sleep(2 ** attempt)
                    continue
                raise APIError("gnomAD returned non-JSON content") from exc
            if not isinstance(data, dict):
                raise APIError("gnomAD returned a non-object JSON response")
            if self._is_transient_graphql_error(data):
                if attempt < 3:
                    self.sleep(2 ** attempt)
                    continue
                # Never persist a transient rate-limit response: otherwise an
                # ordinary or offline rerun would indefinitely replay a failed
                # response from the immutable request cache.
                raise APIError("gnomAD GraphQL rate limit remained active after retries")

            response_sha = hashlib.sha256(response.text.encode()).hexdigest()
            metadata = {
                "url": GNOMAD_API,
                "method": "POST",
                "operation_name": operation_name or "",
                "retrieved_at_utc": utc_now(),
                "request_sha256_utf8": hashlib.sha256(canonical.encode()).hexdigest(),
                "response_sha256_utf8": response_sha,
                "content_type": response.headers.get("Content-Type"),
                "last_modified": response.headers.get("Last-Modified"),
                "etag": response.headers.get("ETag"),
            }
            write_json(payload_path, data)
            write_json(metadata_path, metadata)
            self.manifest[key] = metadata
            return data
        raise APIError("gnomAD request failed after retries")


def _chunks(values: list[str], size: int):
    for start in range(0, len(values), size):
        yield values[start:start + size]


def build_search_query(current_rsids: list[str]) -> tuple[str, dict[str, str]]:
    aliases = {f"r{index:03d}": rsid for index, rsid in enumerate(current_rsids)}
    declarations = ", ".join(f"${alias}: String!" for alias in aliases)
    fields = "\n".join(
        f"  {alias}: variant_search(query: ${alias}, dataset: {GNOMAD_DATASET}) "
        "{ variant_id }"
        for alias in aliases
    )
    return f"query PanelVariantSearch({declarations}) {{\n{fields}\n}}", aliases


DETAIL_FRAGMENT = """
fragment PanelVariantFields on VariantDetails {
  variant_id
  reference_genome
  chrom
  pos
  ref
  alt
  caid
  rsids
  flags
  exome {
    ac an af homozygote_count hemizygote_count filters flags
    faf95 { popmax_population popmax }
    populations { id ac an homozygote_count hemizygote_count }
  }
  genome {
    ac an af homozygote_count hemizygote_count filters flags
    faf95 { popmax_population popmax }
    populations { id ac an homozygote_count hemizygote_count }
  }
  joint {
    ac an homozygote_count hemizygote_count filters
    faf95 { popmax_population popmax }
    populations { id ac an homozygote_count hemizygote_count }
  }
}
""".strip()


def build_detail_query(variant_ids: list[str]) -> tuple[str, dict[str, str]]:
    if len(variant_ids) > DETAIL_CHUNK_SIZE:
        raise ValueError(f"gnomAD detail chunks may contain at most {DETAIL_CHUNK_SIZE} IDs")
    aliases = {f"v{index:03d}": variant_id for index, variant_id in enumerate(variant_ids)}
    declarations = ", ".join(f"${alias}: String!" for alias in aliases)
    fields = "\n".join(
        f"  {alias}: variant(variantId: ${alias}, dataset: {GNOMAD_DATASET}) "
        "{ ...PanelVariantFields }"
        for alias in aliases
    )
    query = (f"query PanelVariantDetails({declarations}) {{\n{fields}\n}}\n"
             f"{DETAIL_FRAGMENT}")
    return query, aliases


def graphql_messages(response: dict, alias: str | None = None) -> list[str]:
    messages = []
    for error in response.get("errors", []) or []:
        if not isinstance(error, dict):
            messages.append(str(error))
            continue
        path = error.get("path") or []
        if alias is None or not path or str(path[0]) == alias:
            messages.append(str(error.get("message", "Unknown GraphQL error")))
    return messages


def _int(value: Any) -> int | None:
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _float(value: Any) -> float | None:
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def calculated_af(ac: int | None, an: int | None) -> float | None:
    return ac / an if ac is not None and an is not None and an > 0 else None


def population_scope(population_id: str) -> str:
    if not population_id:
        return "unlabeled_overall"
    if population_id in {"XX", "XY"}:
        return "overall_sex_karyotype"
    if population_id.endswith("_XX") or population_id.endswith("_XY"):
        return "sex_stratified"
    if population_id.startswith("hgdp:"):
        return "hgdp_subpopulation"
    if population_id.startswith("1kg:"):
        return "one_thousand_genomes_subpopulation"
    return "genetic_ancestry_group"


def normalize_callset(variant_id: str, callset: str,
                      value: Any) -> tuple[dict | None, list[dict], int]:
    if not isinstance(value, dict):
        return None, [], 0
    ac, an = _int(value.get("ac")), _int(value.get("an"))
    faf95 = value.get("faf95") if isinstance(value.get("faf95"), dict) else {}
    total = {
        "variant_id": variant_id,
        "callset": callset,
        "ac": ac,
        "an": an,
        "af": calculated_af(ac, an),
        "reported_af": _float(value.get("af")),
        "homozygote_count": _int(value.get("homozygote_count")),
        "hemizygote_count": _int(value.get("hemizygote_count")),
        "filters": value.get("filters") if isinstance(value.get("filters"), list) else [],
        "flags": value.get("flags") if isinstance(value.get("flags"), list) else [],
        "faf95_popmax_population": faf95.get("popmax_population", ""),
        "faf95_popmax": _float(faf95.get("popmax")),
    }

    population_rows: list[dict] = []
    by_id: dict[str, dict] = {}
    duplicates = 0
    for population in value.get("populations", []) or []:
        if not isinstance(population, dict):
            continue
        population_id = str(population.get("id", ""))
        pop_ac, pop_an = _int(population.get("ac")), _int(population.get("an"))
        row = {
            "variant_id": variant_id,
            "callset": callset,
            "population_id": population_id,
            "population_scope": population_scope(population_id),
            "ac": pop_ac,
            "an": pop_an,
            "af": calculated_af(pop_ac, pop_an),
            "homozygote_count": _int(population.get("homozygote_count")),
            "hemizygote_count": _int(population.get("hemizygote_count")),
        }
        previous = by_id.get(population_id)
        if previous is None:
            by_id[population_id] = row
            population_rows.append(row)
        elif previous == row:
            duplicates += 1
        else:
            raise ValueError(
                f"Conflicting duplicate population {population_id!r} for "
                f"{variant_id} {callset}"
            )
    return total, population_rows, duplicates


def normalize_variant_detail(requested_variant_id: str, detail: Any
                             ) -> tuple[dict, list[dict], list[dict], int]:
    if not isinstance(detail, dict):
        raise ValueError("gnomAD detail is null or not an object")
    returned_id = str(detail.get("variant_id", ""))
    if returned_id != requested_variant_id:
        raise ValueError(
            f"Requested {requested_variant_id} but gnomAD returned {returned_id or 'no ID'}"
        )
    if not VARIANT_ID_RE.fullmatch(returned_id):
        raise ValueError(f"Unexpected gnomAD variant ID format: {returned_id}")
    chrom, pos_text, ref, alt = returned_id.split("-", 3)
    if str(detail.get("reference_genome", "")) != "GRCh38":
        raise ValueError(f"{returned_id} did not return GRCh38 data")
    if (str(detail.get("chrom", "")) != chrom or _int(detail.get("pos")) != int(pos_text)
            or str(detail.get("ref", "")) != ref or str(detail.get("alt", "")) != alt):
        raise ValueError(f"Returned coordinates/alleles do not match {returned_id}")

    identity = {
        "variant_id": returned_id,
        "reference_genome": "GRCh38",
        "chrom": chrom,
        "pos": int(pos_text),
        "ref": ref,
        "alt": alt,
        "caid": detail.get("caid", ""),
        "rsids": list(dict.fromkeys(str(value) for value in (detail.get("rsids") or []))),
        "flags": detail.get("flags") if isinstance(detail.get("flags"), list) else [],
        "source_url": (f"https://gnomad.broadinstitute.org/variant/{returned_id}"
                       f"?dataset={GNOMAD_DATASET}"),
    }
    totals, populations, duplicate_count = [], [], 0
    for callset in ("exome", "genome", "joint"):
        total, rows, duplicates = normalize_callset(returned_id, callset,
                                                    detail.get(callset))
        if total is not None:
            totals.append(total)
        populations.extend(rows)
        duplicate_count += duplicates
    return identity, totals, populations, duplicate_count


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
        "n_search_hits": 0,
        "n_verified_variants": 0,
        "n_detail_missing": 0,
        "n_verification_failed": 0,
        "gnomad_variant_ids": [],
        "error": error,
    }


def _sqlite_value(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, bool):
        return int(value)
    return value


def _columns(fields: list[str], integers: set[str], reals: set[str] | None = None) -> str:
    reals = reals or set()
    return ", ".join(
        f'"{field}" {"INTEGER" if field in integers else "REAL" if field in reals else "TEXT"}'
        for field in fields
    )


def _insert(database: sqlite3.Connection, table: str, fields: list[str],
            rows: list[dict]) -> None:
    quoted = ", ".join(f'"{field}"' for field in fields)
    placeholders = ", ".join("?" for _ in fields)
    database.executemany(
        f'INSERT INTO "{table}" ({quoted}) VALUES ({placeholders})',
        [tuple(_sqlite_value(row.get(field)) for field in fields) for row in rows],
    )


def write_sqlite(path: Path, dbsnp_rows: list[dict], links: list[dict],
                 variants: list[dict], callsets: list[dict], populations: list[dict],
                 statuses: list[dict], manifest: dict, provenance: dict) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    try:
        with database:
            database.execute("PRAGMA foreign_keys = ON")
            database.execute(
                f"CREATE TABLE dbsnp_resolution ({_columns(DBSNP_FIELDS, {'rank', 'was_merged'})}, "
                "PRIMARY KEY (input_rsid))"
            )
            database.execute(
                f"CREATE TABLE variants ({_columns(VARIANT_FIELDS, {'pos'})}, "
                "PRIMARY KEY (variant_id))"
            )
            database.execute(
                f"CREATE TABLE variant_links ({_columns(LINK_FIELDS, {'rank', 'was_merged', 'gnomad_variant_rank', 'detail_available', 'exact_current_rsid_xref'})}, "
                "PRIMARY KEY (rsid, variant_id), FOREIGN KEY (variant_id) REFERENCES variants(variant_id))"
            )
            database.execute(
                f"CREATE TABLE callset_frequencies ({_columns(CALLSET_FIELDS, {'ac', 'an', 'homozygote_count', 'hemizygote_count'}, {'af', 'reported_af', 'faf95_popmax'})}, "
                "PRIMARY KEY (variant_id, callset), FOREIGN KEY (variant_id) REFERENCES variants(variant_id))"
            )
            database.execute(
                f"CREATE TABLE population_frequencies ({_columns(POPULATION_FIELDS, {'ac', 'an', 'homozygote_count', 'hemizygote_count'}, {'af'})}, "
                "PRIMARY KEY (variant_id, callset, population_id), "
                "FOREIGN KEY (variant_id, callset) REFERENCES callset_frequencies(variant_id, callset))"
            )
            status_integers = {
                "rank", "was_merged", "n_search_hits", "n_verified_variants",
                "n_detail_missing", "n_verification_failed",
            }
            database.execute(
                f"CREATE TABLE panel_status ({_columns(STATUS_FIELDS, status_integers)}, "
                "PRIMARY KEY (rsid))"
            )
            database.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT)")
            _insert(database, "dbsnp_resolution", DBSNP_FIELDS, dbsnp_rows)
            _insert(database, "variants", VARIANT_FIELDS, variants)
            _insert(database, "variant_links", LINK_FIELDS, links)
            _insert(database, "callset_frequencies", CALLSET_FIELDS, callsets)
            _insert(database, "population_frequencies", POPULATION_FIELDS, populations)
            _insert(database, "panel_status", STATUS_FIELDS, statuses)
            database.executemany(
                "INSERT INTO metadata VALUES (?, ?)",
                [("manifest", json.dumps(manifest, ensure_ascii=False, sort_keys=True)),
                 ("source_provenance", json.dumps(provenance, ensure_ascii=False,
                                                  sort_keys=True))],
            )
            integrity = database.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise ValueError(f"SQLite integrity check failed: {integrity}")
    finally:
        database.close()
    partial.replace(path)


def write_readme(out: Path, manifest: dict) -> None:
    text = f"""# gnomAD dataset for the QIAGEN 20-rsID panel

This folder is an allele-frequency snapshot from the official gnomAD GraphQL
API for the exact 20 public test rsIDs. dbSNP was used first to resolve merged
identifiers. It contains no patient genotypes.

- Input rsIDs: {manifest['n_input_rsids']}
- Merged input rsIDs resolved by dbSNP: {manifest['n_merged_input_rsids']}
- gnomAD variant links: {manifest['n_verified_links']}
- Unique gnomAD variants: {manifest['n_variants']}
- Input rsIDs with at least one gnomAD variant: {manifest['n_rsids_with_variants']}

## Files

`dbsnp_resolution.tsv` preserves each supplied ID and its current dbSNP ID.
`rsid_gnomad_links.tsv` records the verified rsID-to-gnomAD-allele mapping.
`gnomad_variants.tsv` contains allele identity, while
`callset_frequencies.tsv` keeps exome, genome, and joint totals separate.
`population_frequencies.tsv` expands the population rows and calculates AF as
AC/AN. JSON copies preserve the same tables, `panel_status.tsv` contains all 20
inputs in original order, `raw/` preserves GraphQL and dbSNP responses, and
`source_provenance.json` records retrieval metadata. `gnomad.sqlite` is only
created for a run with no lookup or identity-verification errors.

Population rows overlap: ancestry totals, XX/XY subdivisions, overall XX/XY,
and some HGDP/1000 Genomes rows are not independent groups. Never sum these
rows. For X/Y variants, retain hemizygote counts and do not infer AN as twice
the sample count. Mitochondrial variants use a different gnomAD API schema;
this panel contains none.

An absent gnomAD record is not evidence that an allele is benign or absent in
all humans. Frequencies, filters, flags, coverage, ancestry composition, and
the exact normalized allele must be considered together.

Official sources:
- https://github.com/broadinstitute/gnomad-browser/blob/main/graphql-api/README.md
- https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def run(args: argparse.Namespace, dbsnp_client: Client | None = None,
        graphql_client: GnomadGraphQLClient | None = None) -> int:
    if args.offline and args.refresh:
        raise ValueError("--offline and --refresh are incompatible")
    rsids = select_panel_rsids(args.input, require_exact=True)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # A failed refresh must not leave a previous successful database looking
    # current. Remove only artifacts owned by this generator; user-added files
    # in the dataset directory are left untouched.
    (out / "gnomad.sqlite").unlink(missing_ok=True)
    (out / "gnomad.sqlite.part").unlink(missing_ok=True)
    raw_detail_dir = out / "raw" / "gnomad_variant_details"
    if raw_detail_dir.is_dir():
        for stale_chunk in raw_detail_dir.glob("chunk_*.json"):
            stale_chunk.unlink()
    dbsnp_client = dbsnp_client or Client(args.cache, args.offline, args.refresh)
    graphql_client = graphql_client or GnomadGraphQLClient(
        Path(args.cache) / "gnomad_graphql", args.offline, args.refresh
    )
    limiter = RateLimiter(0 if args.offline else 0.34)
    errors: list[dict] = []

    print(f"Resolving {len(rsids)} supplied IDs with NCBI dbSNP...", flush=True)
    try:
        dbsnp_records, dbsnp_raw = fetch_dbsnp_summaries(dbsnp_client, limiter, rsids)
    except Exception as exc:
        dbsnp_records = {}
        dbsnp_raw = {"requested_rsids": rsids, "error": str(exc)}
    write_json(out / "raw" / "dbsnp_esummary.json", dbsnp_raw)
    dbsnp_rows = [
        dbsnp_resolution_row(rsid, rank, dbsnp_records.get(rsid))
        for rank, rsid in enumerate(rsids, 1)
    ]
    resolution_by_input = {row["input_rsid"]: row for row in dbsnp_rows}
    for row in dbsnp_rows:
        if not row["current_rsid"]:
            errors.append({
                "rank": row["rank"], "rsid": row["input_rsid"],
                "stage": "dbsnp_resolution", "error": row["error"],
            })

    current_rsids = list(dict.fromkeys(
        row["current_rsid"] for row in dbsnp_rows if row["current_rsid"]
    ))
    search_by_current: dict[str, list[str]] = {rsid: [] for rsid in current_rsids}
    search_failures: dict[str, str] = {}
    search_query, search_aliases = build_search_query(current_rsids)
    search_response: dict = {}
    if current_rsids:
        print(f"Searching gnomAD for {len(current_rsids)} current rsIDs...", flush=True)
        try:
            search_response = graphql_client.post(
                search_query, search_aliases, "PanelVariantSearch"
            )
        except Exception as exc:
            search_response = {"data": None, "client_error": str(exc)}
            for current_rsid in current_rsids:
                search_failures[current_rsid] = str(exc)
    write_json(out / "raw" / "gnomad_variant_search.json", {
        "query": search_query, "variables": search_aliases,
        "response": search_response,
    })
    for message in graphql_messages(search_response):
        errors.append({"stage": "gnomad_variant_search", "error": message})
    search_data = search_response.get("data")
    search_data = search_data if isinstance(search_data, dict) else {}
    for alias, current_rsid in search_aliases.items():
        value = search_data.get(alias)
        if not isinstance(value, list):
            message = "; ".join(graphql_messages(search_response, alias)) or \
                search_response.get("client_error") or "Missing gnomAD search result"
            search_failures[current_rsid] = str(message)
            if not any(error.get("stage") == "gnomad_variant_search"
                       and error.get("current_rsid") == current_rsid for error in errors):
                errors.append({"current_rsid": current_rsid,
                               "stage": "gnomad_variant_search", "error": str(message)})
            continue
        candidates = []
        for item in value:
            variant_id = str(item.get("variant_id", "")) if isinstance(item, dict) else ""
            if not VARIANT_ID_RE.fullmatch(variant_id):
                message = f"Invalid gnomAD search variant ID: {variant_id or item!r}"
                errors.append({"current_rsid": current_rsid,
                               "stage": "gnomad_variant_search", "error": message})
                continue
            if variant_id not in candidates:
                candidates.append(variant_id)
        search_by_current[current_rsid] = candidates
        if len(value) >= 100:
            message = "gnomAD variant_search reached its 100-record server cap"
            search_failures[current_rsid] = message
            errors.append({"current_rsid": current_rsid,
                           "stage": "gnomad_variant_search", "error": message})

    unique_variant_ids = list(dict.fromkeys(
        variant_id for current_rsid in current_rsids
        for variant_id in search_by_current[current_rsid]
    ))
    details_by_id: dict[str, dict] = {}
    detail_failures: dict[str, str] = {}
    print(f"Fetching {len(unique_variant_ids)} unique gnomAD variant records...", flush=True)
    for chunk_index, group in enumerate(_chunks(unique_variant_ids, DETAIL_CHUNK_SIZE), 1):
        detail_query, aliases = build_detail_query(group)
        try:
            response = graphql_client.post(detail_query, aliases, "PanelVariantDetails")
        except Exception as exc:
            response = {"data": None, "client_error": str(exc)}
        write_json(out / "raw" / "gnomad_variant_details" /
                   f"chunk_{chunk_index:03d}.json", {
                       "query": detail_query, "variables": aliases, "response": response,
                   })
        for message in graphql_messages(response):
            errors.append({"detail_chunk": chunk_index,
                           "stage": "gnomad_variant_detail", "error": message})
        data = response.get("data")
        data = data if isinstance(data, dict) else {}
        for alias, variant_id in aliases.items():
            detail = data.get(alias)
            if isinstance(detail, dict):
                details_by_id[variant_id] = detail
            else:
                message = "; ".join(graphql_messages(response, alias)) or \
                    response.get("client_error") or "Missing gnomAD variant detail"
                detail_failures[variant_id] = str(message)
                errors.append({"variant_id": variant_id,
                               "stage": "gnomad_variant_detail", "error": str(message)})

    normalized_variants: dict[str, dict] = {}
    totals_by_variant: dict[str, list[dict]] = {}
    populations_by_variant: dict[str, list[dict]] = {}
    duplicate_populations_removed = 0
    for variant_id in unique_variant_ids:
        if variant_id not in details_by_id:
            continue
        try:
            identity, totals, populations, duplicate_count = normalize_variant_detail(
                variant_id, details_by_id[variant_id]
            )
            normalized_variants[variant_id] = identity
            totals_by_variant[variant_id] = totals
            populations_by_variant[variant_id] = populations
            duplicate_populations_removed += duplicate_count
        except Exception as exc:
            detail_failures[variant_id] = str(exc)
            errors.append({"variant_id": variant_id, "stage": "normalize_gnomad_detail",
                           "error": str(exc)})

    links: list[dict] = []
    statuses: list[dict] = []
    verified_variant_ids: set[str] = set()
    for rank, rsid in enumerate(rsids, 1):
        resolution = resolution_by_input[rsid]
        current_rsid = resolution["current_rsid"]
        if not current_rsid:
            statuses.append(empty_status(rank, rsid, resolution, "dbsnp_not_found",
                                         resolution["error"]))
            continue
        if current_rsid in search_failures and not search_by_current[current_rsid]:
            statuses.append(empty_status(rank, rsid, resolution, "search_failed",
                                         search_failures[current_rsid]))
            continue
        candidates = search_by_current[current_rsid]
        if not candidates:
            statuses.append(empty_status(rank, rsid, resolution, "not_in_gnomad"))
            continue

        verified, missing, failed, local_errors = [], 0, 0, []
        for variant_rank, variant_id in enumerate(candidates, 1):
            identity = normalized_variants.get(variant_id)
            detail_available = identity is not None
            exact_xref = bool(identity and current_rsid in identity["rsids"])
            if not detail_available:
                link_status = "detail_missing_or_invalid"
                missing += 1
                local_errors.append(f"{variant_id}: {detail_failures.get(variant_id, 'missing detail')}")
            elif not exact_xref:
                link_status = "rsid_verification_failed"
                failed += 1
                message = (f"{variant_id}: returned rsids {identity['rsids']} do not include "
                           f"{current_rsid}")
                local_errors.append(message)
                errors.append({"rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                               "variant_id": variant_id,
                               "stage": "verify_gnomad_rsid", "error": message})
            else:
                link_status = "verified"
                verified.append(variant_id)
                verified_variant_ids.add(variant_id)
            links.append({
                "rank": rank, "rsid": rsid, "current_rsid": current_rsid,
                "was_merged": current_rsid != rsid,
                "gnomad_variant_rank": variant_rank,
                "variant_id": variant_id,
                "detail_available": detail_available,
                "exact_current_rsid_xref": exact_xref,
                "status": link_status,
                "source_url": (f"https://gnomad.broadinstitute.org/variant/{variant_id}"
                               f"?dataset={GNOMAD_DATASET}"),
            })

        if len(verified) == len(candidates):
            status = "multiple_variants_found" if len(verified) > 1 else "variant_found"
        elif verified:
            status = "partial_variants_found"
        else:
            status = "detail_lookup_or_verification_failed"
        row = empty_status(rank, rsid, resolution, status, "; ".join(local_errors))
        row.update({
            "n_search_hits": len(candidates),
            "n_verified_variants": len(verified),
            "n_detail_missing": missing,
            "n_verification_failed": failed,
            "gnomad_variant_ids": verified,
        })
        statuses.append(row)

    variants = [normalized_variants[variant_id] for variant_id in unique_variant_ids
                if variant_id in verified_variant_ids]
    callsets = [row for variant_id in unique_variant_ids if variant_id in verified_variant_ids
                for row in totals_by_variant.get(variant_id, [])]
    populations = [row for variant_id in unique_variant_ids if variant_id in verified_variant_ids
                   for row in populations_by_variant.get(variant_id, [])]
    verified_links = [row for row in links if row["status"] == "verified"]
    input_path = Path(args.input)
    status_counts = Counter(row["status"] for row in statuses)
    provenance = {
        "ncbi_dbsnp": list(getattr(dbsnp_client, "manifest", {}).values()),
        "gnomad_graphql": list(getattr(graphql_client, "manifest", {}).values()),
    }
    manifest = {
        "finished_at_utc": utc_now(),
        "generator": Path(__file__).name,
        "source_databases": ["NCBI dbSNP", "gnomAD"],
        "gnomad_api": GNOMAD_API,
        "gnomad_dataset": GNOMAD_DATASET,
        "reference_genome": "GRCh38",
        "input": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "n_input_rsids": len(rsids),
        "n_merged_input_rsids": sum(bool(row["was_merged"]) for row in dbsnp_rows),
        "n_search_hits": sum(len(values) for values in search_by_current.values()),
        "n_verified_links": len(verified_links),
        "n_variants": len(variants),
        "n_rsids_with_variants": sum(row["n_verified_variants"] > 0 for row in statuses),
        "n_callset_rows": len(callsets),
        "n_population_rows": len(populations),
        "n_duplicate_population_rows_removed": duplicate_populations_removed,
        "n_lookup_or_verification_errors": len(errors),
        "panel_status_counts": dict(sorted(status_counts.items())),
        "sqlite_created": False,
        "notes": [
            "Supplied rsIDs remain primary panel keys; merged/current dbSNP IDs are separate",
            "variant_search discovers all gnomAD-present alleles and gnomAD-normalized indels",
            "Every detail record must cross-reference the resolved current rsID",
            "Exome, genome, and joint callsets remain separate",
            "AF is calculated as AC/AN; population strata overlap and must not be summed",
            "Absence from gnomAD is not evidence of benignity",
            "The input contains public loci, not an individual's genotype",
        ],
    }

    write_json(out / "dbsnp_resolution.json", dbsnp_rows)
    write_tsv(out / "dbsnp_resolution.tsv", dbsnp_rows, DBSNP_FIELDS)
    write_json(out / "rsid_gnomad_links.json", links)
    write_tsv(out / "rsid_gnomad_links.tsv", links, LINK_FIELDS)
    write_json(out / "gnomad_variants.json", variants)
    write_tsv(out / "gnomad_variants.tsv", variants, VARIANT_FIELDS)
    write_json(out / "callset_frequencies.json", callsets)
    write_tsv(out / "callset_frequencies.tsv", callsets, CALLSET_FIELDS)
    write_json(out / "population_frequencies.json", populations)
    write_tsv(out / "population_frequencies.tsv", populations, POPULATION_FIELDS)
    write_json(out / "panel_status.json", statuses)
    write_tsv(out / "panel_status.tsv", statuses, STATUS_FIELDS)
    write_json(out / "source_provenance.json", provenance)
    write_json(out / "errors.json", errors)
    write_json(out / "run_manifest.json", manifest)
    write_readme(out, manifest)
    if not errors:
        sqlite_manifest = {**manifest, "sqlite_created": True}
        write_sqlite(out / "gnomad.sqlite", dbsnp_rows, verified_links, variants,
                     callsets, populations, statuses, sqlite_manifest, provenance)
        manifest = sqlite_manifest
        write_json(out / "run_manifest.json", manifest)

    print(f"gnomAD dataset complete: {len(rsids)} rsIDs, {len(verified_links)} links, "
          f"{len(variants)} unique variants, {len(errors)} errors. See {out}")
    return 1 if errors else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="Exact 20-rsID QIAGEN TSV in the original order")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Dataset output directory")
    parser.add_argument("--cache", default=".cache", help="Shared public-API response cache")
    parser.add_argument("--offline", action="store_true", help="Use cache only")
    parser.add_argument("--refresh", action="store_true", help="Refresh cached requests")
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
