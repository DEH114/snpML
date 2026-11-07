"""Shared I/O, conservative SNV normalization and cached HTTP for BioHack.
All genomic coordinates are 1-based, forward-strand GRCh38.
No credentials are required for the public endpoints used here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlsplit
import requests

ENSEMBL = "https://rest.ensembl.org"
DNA = set("ACGT")
CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: str | Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def write_tsv(path: str | Path, rows: list[dict], fields: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})


class APIError(RuntimeError):
    def __init__(self, message: str, status: int | None = None):
        super().__init__(message)
        self.status = status


class Client:
    """GET-only, serial requests; caches body + URL + retrieval date + SHA256.
    HTTP 429/5xx and connection errors are retried. Other errors fail explicitly.
    --offline never makes network calls. Cache is immutable unless --refresh is used.
    """
    def __init__(self, cache: str | Path = ".cache", offline: bool = False,
                 refresh: bool = False, timeout: int = 45):
        if offline and refresh:
            raise ValueError("--offline and --refresh cannot be used together")
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline, self.refresh, self.timeout = offline, refresh, timeout
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "BioHack-SNP-Toolkit/1.0 (research prototype)"
        self.manifest: dict[str, dict] = {}
        self.last_request: dict[str, float] = {}

    def text(self, url: str, params: dict | None = None, accept: str = "application/json") -> str:
        prepared = requests.Request("GET", url, params=params).prepare().url
        key = hashlib.sha256((prepared + "\n" + accept).encode()).hexdigest()
        payload, meta_path = self.cache / f"{key}.txt", self.cache / f"{key}.json"
        if payload.exists() and meta_path.exists() and not self.refresh:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            self.manifest[key] = meta
            return payload.read_text(encoding="utf-8")
        if self.offline:
            raise APIError(f"Offline cache miss: {prepared}")
        host = urlsplit(prepared).netloc
        for attempt in range(4):
            delay = 0.25 - (time.monotonic() - self.last_request.get(host, 0))
            if delay > 0:
                time.sleep(delay)
            self.last_request[host] = time.monotonic()
            try:
                response = self.session.get(prepared, headers={"Accept": accept},
                                            timeout=(10, self.timeout))
            except requests.RequestException as exc:
                if attempt == 3:
                    raise APIError(f"Network error for {prepared}: {exc}") from exc
                time.sleep(2 ** attempt)
                continue
            if response.status_code in (429, 500, 502, 503, 504):
                if attempt < 3:
                    try:
                        wait = min(float(response.headers.get("Retry-After", 2 ** attempt)), 60)
                    except ValueError:
                        wait = 2 ** attempt
                    time.sleep(wait)
                    continue
            if not response.ok:
                raise APIError(f"HTTP {response.status_code}: {prepared}: {response.text[:350]}",
                               response.status_code)
            text = response.text
            if accept == "application/json":
                try:
                    json.loads(text)
                except json.JSONDecodeError as exc:
                    raise APIError(f"Expected JSON but received another format: {prepared}") from exc
            meta = {"url": response.url, "retrieved_at_utc": utc_now(),
                    "sha256_utf8": hashlib.sha256(text.encode()).hexdigest(),
                    "content_type": response.headers.get("Content-Type"),
                    "last_modified": response.headers.get("Last-Modified"),
                    "etag": response.headers.get("ETag")}
            payload.write_text(text, encoding="utf-8")
            write_json(meta_path, meta)
            self.manifest[key] = meta
            return text
        raise APIError(f"Request failed: {prepared}")

    def json(self, url: str, params: dict | None = None) -> Any:
        return json.loads(self.text(url, params=params))

    def ensembl(self, path: str, params: dict | None = None) -> Any:
        p = {"content-type": "application/json", **(params or {})}
        return self.json(ENSEMBL + path, p)


def revcomp(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str = ""
    assembly: str = "GRCh38"

    def __post_init__(self):
        if self.assembly != "GRCh38":
            raise ValueError("Only GRCh38 is supported; do not relabel GRCh37 coordinates")
        if self.chrom not in CHROMS or self.pos < 1:
            raise ValueError("Require a primary chromosome and positive 1-based position")
        if self.ref not in DNA or self.alt not in DNA or self.ref == self.alt:
            raise ValueError("Only single-base A/C/G/T substitutions with REF != ALT are supported")

    @property
    def key(self) -> str:
        return f"{self.assembly}:{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    def row(self) -> dict:
        return {"variant_key": self.key, **asdict(self)}


def read_inputs(path: str | Path) -> list[dict[str, str]]:
    """TSV or one rsID per line. TSV accepts rsid[,alt] OR chrom,pos,ref,alt[,assembly]."""
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = [x for x in text.splitlines() if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        raise ValueError("Input file is empty")
    if all(re.fullmatch(r"rs\d+", x.strip()) for x in lines):
        return [{"rsid": x.strip()} for x in lines]
    reader = csv.DictReader(lines, delimiter="\t")
    if not reader.fieldnames or not ({"rsid"} <= set(reader.fieldnames) or
                                    {"chrom", "pos", "ref", "alt"} <= set(reader.fieldnames)):
        raise ValueError("Use TSV headers: rsid,alt OR chrom,pos,ref,alt,assembly")
    return [{k: (v or "").strip() for k, v in row.items() if k is not None} for row in reader]


def sequence(client: Client, chrom: str, start: int, end: int) -> str:
    data = client.ensembl(f"/sequence/region/human/{chrom}:{start}..{end}:1",
                          {"coord_system_version": "GRCh38"})
    seq = data.get("seq", "").upper()
    if len(seq) != end - start + 1:
        raise ValueError(f"Unexpected reference sequence length at {chrom}:{start}-{end}")
    return seq


def mapping_candidates(data: dict, rsid: str, alt: str = "") -> list[Variant]:
    candidates: dict[tuple, set[str]] = {}
    for m in data.get("mappings", []):
        if m.get("assembly_name") != "GRCh38" or str(m.get("seq_region_name")) not in CHROMS:
            continue
        if m.get("start") != m.get("end"):
            continue
        alleles = str(m.get("allele_string", "")).upper().split("/")
        if int(m.get("strand", 1)) == -1:
            alleles = [revcomp(a) for a in alleles]
        if not alleles or alleles[0] not in DNA:
            continue
        key = (str(m["seq_region_name"]), int(m["start"]), alleles[0])
        candidates.setdefault(key, set()).update(a for a in alleles[1:] if a in DNA and a != alleles[0])
    candidates = {k: v for k, v in candidates.items() if v}
    if len(candidates) != 1:
        raise ValueError(f"{rsid}: expected one unambiguous GRCh38 primary SNV mapping, got {len(candidates)}")
    (chrom, pos, ref), alts = next(iter(candidates.items()))
    if alt:
        if alt not in alts:
            raise ValueError(f"{rsid}: ALT {alt} not in forward-strand SNV alleles {sorted(alts)}")
        alts = {alt}
    return [Variant(chrom, pos, ref, a, rsid) for a in sorted(alts)]


def resolve(client: Client, row: dict) -> list[Variant]:
    assembly = row.get("assembly") or "GRCh38"
    if assembly != "GRCh38":
        raise ValueError("This implementation accepts GRCh38 only")
    coordinates = [row.get(k, "") for k in ("chrom", "pos", "ref")]
    if all(coordinates):
        chrom = re.sub(r"^chr", "", row["chrom"], flags=re.I).upper()
        chrom = "MT" if chrom == "M" else chrom
        variants = [Variant(chrom, int(row["pos"]), row["ref"].upper(),
                            row.get("alt", "").upper(), row.get("rsid", ""), assembly)]
    elif any(coordinates):
        raise ValueError("Partially specified coordinates; provide chrom,pos,ref,alt together")
    else:
        rsid = row.get("rsid", "")
        if not re.fullmatch(r"rs\d+", rsid):
            raise ValueError(f"Not an rsID: {rsid!r}")
        data = client.ensembl(f"/variation/human/{rsid}")
        variants = mapping_candidates(data, rsid, row.get("alt", "").upper())
    for v in variants:
        actual = sequence(client, v.chrom, v.pos, v.pos)
        if actual != v.ref:
            raise ValueError(f"REF mismatch for {v.key}: genome has {actual}; check build/strand")
    return variants


def vep(client: Client, v: Variant, alphamissense: bool = True) -> list[dict]:
    params = {"mane": 1, "canonical": 1, "protein": 1, "uniprot": 1,
              "hgvs": 1, "domains": 1}
    if alphamissense:
        params["AlphaMissense"] = 1
    data = client.ensembl(f"/vep/human/region/{v.chrom}:{v.pos}-{v.pos}:1/{v.alt}", params)
    if not isinstance(data, list) or not data:
        raise ValueError("VEP returned no annotation")
    result = []
    for record in data:
        if (record.get("assembly_name") != "GRCh38" or
            str(record.get("seq_region_name")) != v.chrom or record.get("start") != v.pos):
            continue
        alleles = str(record.get("allele_string", "")).split("/")
        if alleles and alleles[0] == v.ref and v.alt in alleles[1:]:
            result.append(record)
    if not result:
        raise ValueError(f"VEP output did not match exact genomic allele: {v.key}")
    return result


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", required=True, help="TSV input or one rsID per line")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--cache", default=".cache", help="Shared response-cache directory")
    parser.add_argument("--offline", action="store_true", help="Use cache only, never access network")
    parser.add_argument("--refresh", action="store_true", help="Ignore old cached responses")


def save_run(out: Path, client: Client, args: argparse.Namespace, errors: list, extra: dict) -> None:
    write_json(out / "errors.json", errors)
    write_json(out / "run_manifest.json", {"finished_at_utc": utc_now(), "assembly": "GRCh38",
        "arguments": vars(args),
        "input_sha256": hashlib.sha256(Path(args.input).read_bytes()).hexdigest(),
        "sources": list(client.manifest.values()), **extra})
