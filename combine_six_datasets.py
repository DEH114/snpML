"""Combine the six normalized rsID datasets into one auditable dataset.

The combiner is deliberately offline: it reads only the named normalized JSON
artifacts and manifests produced by the six generators.  It never scans raw
response directories.  Every output record represents exactly one input rsID;
one-to-many data remain nested under their source instead of being flattened
through a Cartesian join.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from challenge_batch import DEFAULT_INPUT, select_panel_rsids
from common import utc_now, write_json, write_tsv


ROOT = Path(__file__).resolve().parent
DEFAULT_DATASET_ROOT = ROOT / "results" / "qiagen_20_snp_dataset"
SOURCE_NAMES = (
    "alpha_fold", "jaspar", "clinvar", "gwas_catalog", "gnomad", "gtr",
)
SOURCE_DIRS = {name: f"{name}_dataset" for name in SOURCE_NAMES}
OUTPUT_FILENAMES = (
    "combined_data.json",
    "combined_summary.tsv",
    "combined_data.sqlite",
    "combined_data.sqlite.part",
    "run_manifest.json",
    "README.md",
)

SUMMARY_FIELDS = [
    "rank", "rsid", "current_rsid", "was_merged",
    "alpha_fold_status", "n_alpha_fold_alleles", "n_protein_annotations",
    "jaspar_status", "n_jaspar_alleles", "n_motif_hits",
    "clinvar_status", "n_clinvar_links", "n_clinvar_direct_records",
    "n_clinvar_related_set_records", "gwas_catalog_status",
    "n_gwas_associations", "gnomad_status", "n_gnomad_variant_links",
    "n_gnomad_verified_variants", "gtr_status", "n_gtr_test_matches",
    "alpha_fold_json", "jaspar_json", "clinvar_json",
    "gwas_catalog_json", "gnomad_json", "gtr_json",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class SourceReader:
    """Read explicit source files and retain content hashes for provenance."""

    def __init__(self, source: str, directory: Path, dataset_root: Path):
        self.source = source
        self.directory = directory
        self.dataset_root = dataset_root
        self.files: list[dict[str, Any]] = []

    def load(self, filename: str, expected: type | tuple[type, ...]) -> Any:
        path = self.directory / filename
        if not path.is_file():
            raise FileNotFoundError(f"{self.source}: required normalized file missing: {path}")
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"{self.source}: cannot read valid JSON from {path}: {exc}") from exc
        if not isinstance(value, expected):
            names = (expected.__name__ if isinstance(expected, type)
                     else "/".join(item.__name__ for item in expected))
            raise ValueError(f"{self.source}: {filename} must contain {names}")
        self.files.append({
            "relative_path": self._relative(path),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
        return value

    def load_optional(self, filename: str, expected: type | tuple[type, ...]) -> Any | None:
        path = self.directory / filename
        if not path.exists():
            return None
        return self.load(filename, expected)

    def _relative(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.dataset_root.resolve()).as_posix()
        except ValueError:
            return str(path.resolve())


def list_of_dicts(value: Any, label: str) -> list[dict]:
    if not isinstance(value, list) or any(not isinstance(row, dict) for row in value):
        raise ValueError(f"{label} must be a JSON list of objects")
    return value


def validate_panel(rows: list[dict], rsids: list[str], key: str, label: str) -> None:
    observed = [str(row.get(key, "")) for row in rows]
    if observed != rsids:
        raise ValueError(
            f"{label} does not cover the exact input rsIDs in order: "
            f"expected={rsids}, observed={observed}"
        )


def validate_relation_rsids(rows: Iterable[dict], rsids: list[str], label: str,
                            key: str = "rsid") -> None:
    allowed = set(rsids)
    invalid = sorted({str(row.get(key, "")) for row in rows
                      if str(row.get(key, "")) not in allowed})
    if invalid:
        raise ValueError(f"{label} contains rsIDs outside the input panel: {invalid}")


def unique_index(rows: Iterable[dict], key: str, label: str) -> dict[str, dict]:
    result: dict[str, dict] = {}
    for row in rows:
        value = str(row.get(key, ""))
        if not value:
            raise ValueError(f"{label} contains a row without {key}")
        if value in result:
            raise ValueError(f"{label} contains duplicate {key}={value}")
        result[value] = row
    return result


def group_by(rows: Iterable[dict], key: str) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get(key, ""))].append(row)
    return grouped


def panel_index(rows: list[dict], key: str = "rsid") -> dict[str, dict]:
    return {str(row[key]): row for row in rows}


def source_provenance(reader: SourceReader, manifests: dict[str, Any],
                      errors: list[dict]) -> dict:
    return {
        "directory": reader._relative(reader.directory),
        "manifests": manifests,
        "reported_errors": errors,
        "input_files": reader.files,
    }


def load_alpha_fold(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("alpha_fold", directory, root)
    manifest = reader.load("dataset_manifest.json", dict)
    branch_manifest = reader.load_optional("run_manifest.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "alpha_fold/errors.json")
    summaries = list_of_dicts(reader.load("dataset_summary.json", list),
                              "alpha_fold/dataset_summary.json")
    resolutions = list_of_dicts(reader.load("resolution_details.json", list),
                                "alpha_fold/resolution_details.json")
    variants = list_of_dicts(reader.load("variant_status.json", list),
                             "alpha_fold/variant_status.json")
    proteins = list_of_dicts(reader.load("protein_results.json", list),
                             "alpha_fold/protein_results.json")
    validate_panel(summaries, rsids, "rsid", "alpha_fold dataset summary")
    validate_panel(resolutions, rsids, "rsid", "alpha_fold resolution details")
    validate_relation_rsids(variants, rsids, "alpha_fold variant status")
    validate_relation_rsids(proteins, rsids, "alpha_fold protein results")
    variant_index = unique_index(variants, "variant_key", "alpha_fold variant status")
    proteins_by_variant = group_by(proteins, "variant_key")
    unknown = sorted(set(proteins_by_variant) - set(variant_index))
    if unknown:
        raise ValueError(f"alpha_fold protein results reference unknown variants: {unknown}")
    variants_by_rsid = group_by(variants, "rsid")
    summary_by_rsid = panel_index(summaries)
    resolution_by_rsid = panel_index(resolutions)
    records = {}
    for rsid in rsids:
        alleles = []
        for variant in variants_by_rsid.get(rsid, []):
            allele = copy.deepcopy(variant)
            allele["protein_annotations"] = copy.deepcopy(
                proteins_by_variant.get(str(variant["variant_key"]), [])
            )
            alleles.append(allele)
        records[rsid] = {
            "summary": copy.deepcopy(summary_by_rsid[rsid]),
            "resolution": copy.deepcopy(resolution_by_rsid[rsid]),
            "alleles": alleles,
        }
    provenance = source_provenance(
        reader,
        {"dataset_manifest": manifest, "branch_run_manifest": branch_manifest},
        errors,
    )
    return records, provenance


def load_jaspar(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("jaspar", directory, root)
    manifest = reader.load("dataset_manifest.json", dict)
    branch_manifest = reader.load_optional("run_manifest.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "jaspar/errors.json")
    summaries = list_of_dicts(reader.load("dataset_summary.json", list),
                              "jaspar/dataset_summary.json")
    resolutions = list_of_dicts(reader.load("resolution_details.json", list),
                                "jaspar/resolution_details.json")
    variants = list_of_dicts(reader.load("variant_status.json", list),
                             "jaspar/variant_status.json")
    motifs = list_of_dicts(reader.load("motif_results.json", list),
                           "jaspar/motif_results.json")
    profiles = list_of_dicts(reader.load("motif_profiles.json", list),
                             "jaspar/motif_profiles.json")
    validate_panel(summaries, rsids, "rsid", "jaspar dataset summary")
    validate_panel(resolutions, rsids, "rsid", "jaspar resolution details")
    validate_relation_rsids(variants, rsids, "jaspar variant status")
    validate_relation_rsids(motifs, rsids, "jaspar motif results")
    variant_index = unique_index(variants, "variant_key", "jaspar variant status")
    profile_index = unique_index(profiles, "matrix_id", "jaspar motif profiles")
    motifs_by_variant = group_by(motifs, "variant_key")
    unknown_variants = sorted(set(motifs_by_variant) - set(variant_index))
    if unknown_variants:
        raise ValueError(f"jaspar motif results reference unknown variants: {unknown_variants}")
    unknown_profiles = sorted({str(row.get("matrix_id", "")) for row in motifs}
                              - set(profile_index))
    if unknown_profiles:
        raise ValueError(f"jaspar motif results reference unknown profiles: {unknown_profiles}")
    variants_by_rsid = group_by(variants, "rsid")
    summary_by_rsid = panel_index(summaries)
    resolution_by_rsid = panel_index(resolutions)
    records = {}
    for rsid in rsids:
        alleles, used_profiles = [], []
        seen_profiles: set[str] = set()
        for variant in variants_by_rsid.get(rsid, []):
            allele = copy.deepcopy(variant)
            hits = copy.deepcopy(motifs_by_variant.get(str(variant["variant_key"]), []))
            allele["motif_hits"] = hits
            alleles.append(allele)
            for hit in hits:
                matrix_id = str(hit["matrix_id"])
                if matrix_id not in seen_profiles:
                    seen_profiles.add(matrix_id)
                    used_profiles.append(copy.deepcopy(profile_index[matrix_id]))
        records[rsid] = {
            "summary": copy.deepcopy(summary_by_rsid[rsid]),
            "resolution": copy.deepcopy(resolution_by_rsid[rsid]),
            "alleles": alleles,
            "referenced_motif_profiles": used_profiles,
        }
    provenance = source_provenance(
        reader,
        {"dataset_manifest": manifest, "branch_run_manifest": branch_manifest},
        errors,
    )
    return records, provenance


def load_clinvar(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("clinvar", directory, root)
    manifest = reader.load("run_manifest.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "clinvar/errors.json")
    resolutions = list_of_dicts(reader.load("dbsnp_resolution.json", list),
                                "clinvar/dbsnp_resolution.json")
    statuses = list_of_dicts(reader.load("panel_status.json", list),
                             "clinvar/panel_status.json")
    links = list_of_dicts(reader.load("rsid_clinvar_links.json", list),
                          "clinvar/rsid_clinvar_links.json")
    direct = list_of_dicts(reader.load("clinvar_records.json", list),
                           "clinvar/clinvar_records.json")
    related = list_of_dicts(reader.load("related_set_records.json", list),
                            "clinvar/related_set_records.json")
    traits = list_of_dicts(reader.load("traits.json", list), "clinvar/traits.json")
    validate_panel(resolutions, rsids, "input_rsid", "clinvar dbSNP resolution")
    validate_panel(statuses, rsids, "rsid", "clinvar panel status")
    for label, rows in (("links", links), ("direct records", direct),
                        ("related records", related), ("traits", traits)):
        validate_relation_rsids(rows, rsids, f"clinvar {label}")
    known_records = {(str(row.get("rsid", "")), str(row.get("clinvar_uid", "")),
                      str(row.get("record_scope", ""))) for row in direct + related}
    orphan_traits = sorted({(str(row.get("rsid", "")), str(row.get("clinvar_uid", "")),
                             str(row.get("record_scope", ""))) for row in traits}
                           - known_records)
    if orphan_traits:
        raise ValueError(f"clinvar traits reference unknown records: {orphan_traits}")
    traits_by_record = group_by(
        ({**row, "_join_key": "\t".join((str(row.get("rsid", "")),
                                           str(row.get("clinvar_uid", "")),
                                           str(row.get("record_scope", ""))))}
         for row in traits), "_join_key"
    )
    def records_with_traits(rows: list[dict]) -> dict[str, list[dict]]:
        grouped: dict[str, list[dict]] = defaultdict(list)
        for row in rows:
            item = copy.deepcopy(row)
            key = "\t".join((str(row.get("rsid", "")), str(row.get("clinvar_uid", "")),
                              str(row.get("record_scope", ""))))
            item["traits"] = [{k: v for k, v in trait.items() if k != "_join_key"}
                              for trait in copy.deepcopy(traits_by_record.get(key, []))]
            grouped[str(row["rsid"])].append(item)
        return grouped
    direct_by_rsid = records_with_traits(direct)
    related_by_rsid = records_with_traits(related)
    status_by_rsid = panel_index(statuses)
    resolution_by_rsid = panel_index(resolutions, "input_rsid")
    links_by_rsid = group_by(links, "rsid")
    records = {rsid: {
        "status": copy.deepcopy(status_by_rsid[rsid]),
        "dbsnp_resolution": copy.deepcopy(resolution_by_rsid[rsid]),
        "links": copy.deepcopy(links_by_rsid.get(rsid, [])),
        "direct_records": direct_by_rsid.get(rsid, []),
        "related_set_records": related_by_rsid.get(rsid, []),
    } for rsid in rsids}
    return records, source_provenance(reader, {"run_manifest": manifest}, errors)


def load_gwas(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("gwas_catalog", directory, root)
    manifest = reader.load("run_manifest.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "gwas_catalog/errors.json")
    resolutions = list_of_dicts(reader.load("dbsnp_resolution.json", list),
                                "gwas_catalog/dbsnp_resolution.json")
    statuses = list_of_dicts(reader.load("panel_status.json", list),
                             "gwas_catalog/panel_status.json")
    associations = list_of_dicts(reader.load("rsid_associations.json", list),
                                 "gwas_catalog/rsid_associations.json")
    source_records = list_of_dicts(reader.load("matched_source_records.json", list),
                                   "gwas_catalog/matched_source_records.json")
    validate_panel(resolutions, rsids, "input_rsid", "gwas_catalog dbSNP resolution")
    validate_panel(statuses, rsids, "rsid", "gwas_catalog panel status")
    validate_relation_rsids(associations, rsids, "gwas_catalog associations")
    source_index = unique_index(source_records, "source_row_number",
                                "gwas_catalog matched source records")
    unknown = sorted({str(row.get("source_row_number", "")) for row in associations}
                     - set(source_index))
    if unknown:
        raise ValueError(f"gwas_catalog associations reference unknown source records: {unknown}")
    associations_by_rsid: dict[str, list[dict]] = defaultdict(list)
    for association in associations:
        item = copy.deepcopy(association)
        source_record = source_index[str(association["source_row_number"])]
        if (str(source_record.get("source_record_sha256", ""))
                != str(association.get("source_record_sha256", ""))):
            raise ValueError(
                "gwas_catalog association/source digest mismatch at source row "
                f"{association['source_row_number']}"
            )
        item["source_record"] = copy.deepcopy(source_record)
        associations_by_rsid[str(association["rsid"])].append(item)
    status_by_rsid = panel_index(statuses)
    resolution_by_rsid = panel_index(resolutions, "input_rsid")
    records = {rsid: {
        "status": copy.deepcopy(status_by_rsid[rsid]),
        "dbsnp_resolution": copy.deepcopy(resolution_by_rsid[rsid]),
        "associations": associations_by_rsid.get(rsid, []),
    } for rsid in rsids}
    return records, source_provenance(reader, {"run_manifest": manifest}, errors)


def load_gnomad(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("gnomad", directory, root)
    manifest = reader.load("run_manifest.json", dict)
    provenance_manifest = reader.load("source_provenance.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "gnomad/errors.json")
    resolutions = list_of_dicts(reader.load("dbsnp_resolution.json", list),
                                "gnomad/dbsnp_resolution.json")
    statuses = list_of_dicts(reader.load("panel_status.json", list),
                             "gnomad/panel_status.json")
    links = list_of_dicts(reader.load("rsid_gnomad_links.json", list),
                          "gnomad/rsid_gnomad_links.json")
    variants = list_of_dicts(reader.load("gnomad_variants.json", list),
                             "gnomad/gnomad_variants.json")
    callsets = list_of_dicts(reader.load("callset_frequencies.json", list),
                             "gnomad/callset_frequencies.json")
    populations = list_of_dicts(reader.load("population_frequencies.json", list),
                                "gnomad/population_frequencies.json")
    validate_panel(resolutions, rsids, "input_rsid", "gnomad dbSNP resolution")
    validate_panel(statuses, rsids, "rsid", "gnomad panel status")
    validate_relation_rsids(links, rsids, "gnomad links")
    variant_index = unique_index(variants, "variant_id", "gnomad variants")
    linked_ids = {str(row.get("variant_id", "")) for row in links}
    orphan_variants = sorted(set(variant_index) - linked_ids)
    if orphan_variants:
        raise ValueError(f"gnomad variants have no rsID link: {orphan_variants}")
    for link in links:
        if link.get("detail_available") and str(link.get("variant_id", "")) not in variant_index:
            raise ValueError(f"gnomad detailed link has no variant: {link.get('variant_id')}")
    callsets_by_variant = group_by(callsets, "variant_id")
    unknown_callsets = sorted(set(callsets_by_variant) - set(variant_index))
    if unknown_callsets:
        raise ValueError(f"gnomad callsets reference unknown variants: {unknown_callsets}")
    population_by_callset: dict[tuple[str, str], list[dict]] = defaultdict(list)
    callset_keys = {(str(row.get("variant_id", "")), str(row.get("callset", "")))
                    for row in callsets}
    for row in populations:
        key = (str(row.get("variant_id", "")), str(row.get("callset", "")))
        if key not in callset_keys:
            raise ValueError(f"gnomad population frequency references unknown callset: {key}")
        population_by_callset[key].append(row)
    enriched_variants: dict[str, dict] = {}
    for variant_id, variant in variant_index.items():
        item = copy.deepcopy(variant)
        item["callsets"] = []
        for callset in callsets_by_variant.get(variant_id, []):
            enriched = copy.deepcopy(callset)
            enriched["populations"] = copy.deepcopy(population_by_callset.get(
                (variant_id, str(callset.get("callset", ""))), []
            ))
            item["callsets"].append(enriched)
        enriched_variants[variant_id] = item
    links_by_rsid = group_by(links, "rsid")
    status_by_rsid = panel_index(statuses)
    resolution_by_rsid = panel_index(resolutions, "input_rsid")
    records = {}
    for rsid in rsids:
        entries = []
        for link in links_by_rsid.get(rsid, []):
            entry = copy.deepcopy(link)
            entry["variant"] = copy.deepcopy(enriched_variants.get(str(link["variant_id"])))
            entries.append(entry)
        records[rsid] = {
            "status": copy.deepcopy(status_by_rsid[rsid]),
            "dbsnp_resolution": copy.deepcopy(resolution_by_rsid[rsid]),
            "variant_links": entries,
        }
    return records, source_provenance(
        reader, {"run_manifest": manifest, "source_provenance": provenance_manifest}, errors
    )


def load_gtr(directory: Path, root: Path, rsids: list[str]) -> tuple[dict, dict]:
    reader = SourceReader("gtr", directory, root)
    manifest = reader.load("run_manifest.json", dict)
    errors = list_of_dicts(reader.load("errors.json", list), "gtr/errors.json")
    resolutions = list_of_dicts(reader.load("dbsnp_resolution.json", list),
                                "gtr/dbsnp_resolution.json")
    statuses = list_of_dicts(reader.load("panel_status.json", list),
                             "gtr/panel_status.json")
    matches = list_of_dicts(reader.load("gtr_test_matches.json", list),
                            "gtr/gtr_test_matches.json")
    tests = list_of_dicts(reader.load("gtr_tests.json", list), "gtr/gtr_tests.json")
    conditions = list_of_dicts(reader.load("conditions.json", list), "gtr/conditions.json")
    targets = list_of_dicts(reader.load("variant_targets.json", list),
                            "gtr/variant_targets.json")
    validate_panel(resolutions, rsids, "input_rsid", "gtr dbSNP resolution")
    validate_panel(statuses, rsids, "rsid", "gtr panel status")
    validate_relation_rsids(matches, rsids, "gtr test matches")
    test_index = unique_index(tests, "gtr_uid", "gtr tests")
    matched_ids = {str(row.get("gtr_uid", "")) for row in matches}
    orphan_tests = sorted(set(test_index) - matched_ids)
    if orphan_tests:
        raise ValueError(f"gtr tests have no rsID match: {orphan_tests}")
    for match in matches:
        if str(match.get("gtr_uid", "")) not in test_index:
            raise ValueError(f"gtr match references unknown test: {match.get('gtr_uid')}")
    conditions_by_test = group_by(conditions, "gtr_uid")
    targets_by_test = group_by(targets, "gtr_uid")
    unknown_children = sorted((set(conditions_by_test) | set(targets_by_test)) - set(test_index))
    if unknown_children:
        raise ValueError(f"gtr child records reference unknown tests: {unknown_children}")
    matches_by_rsid = group_by(matches, "rsid")
    status_by_rsid = panel_index(statuses)
    resolution_by_rsid = panel_index(resolutions, "input_rsid")
    records = {}
    for rsid in rsids:
        entries = []
        for match in matches_by_rsid.get(rsid, []):
            entry = copy.deepcopy(match)
            uid = str(match["gtr_uid"])
            test = copy.deepcopy(test_index[uid])
            # Preserve the test summary's original ``conditions`` value and
            # attach the separately normalized condition relation under a
            # distinct name.
            test["condition_records"] = copy.deepcopy(conditions_by_test.get(uid, []))
            test["variant_targets"] = copy.deepcopy(targets_by_test.get(uid, []))
            entry["test"] = test
            entries.append(entry)
        records[rsid] = {
            "status": copy.deepcopy(status_by_rsid[rsid]),
            "dbsnp_resolution": copy.deepcopy(resolution_by_rsid[rsid]),
            "test_matches": entries,
        }
    return records, source_provenance(reader, {"run_manifest": manifest}, errors)


LOADERS = {
    "alpha_fold": load_alpha_fold,
    "jaspar": load_jaspar,
    "clinvar": load_clinvar,
    "gwas_catalog": load_gwas,
    "gnomad": load_gnomad,
    "gtr": load_gtr,
}


def status_value(source: dict | None) -> str:
    if source is None:
        return "source_not_available"
    summary = source.get("summary") or source.get("status") or {}
    if "status" in summary:
        return str(summary.get("status", ""))
    if "protein_status" in summary and summary.get("protein_status") != "not_run":
        return str(summary.get("protein_status", ""))
    if "jaspar_status" in summary:
        return str(summary.get("jaspar_status", ""))
    return ""


def count_nested(source: dict | None, field: str) -> int:
    return len(source.get(field, [])) if source else 0


def build_summary(records: list[dict]) -> list[dict]:
    rows = []
    for record in records:
        alpha = record.get("alpha_fold")
        jaspar = record.get("jaspar")
        clinvar = record.get("clinvar")
        gwas = record.get("gwas_catalog")
        gnomad = record.get("gnomad")
        gtr = record.get("gtr")
        alpha_alleles = alpha.get("alleles", []) if alpha else []
        jaspar_alleles = jaspar.get("alleles", []) if jaspar else []
        gnomad_links = gnomad.get("variant_links", []) if gnomad else []
        rows.append({
            "rank": record["rank"], "rsid": record["rsid"],
            "current_rsid": record["identity"]["current_rsid"],
            "was_merged": record["identity"]["was_merged"],
            "alpha_fold_status": status_value(alpha),
            "n_alpha_fold_alleles": len(alpha_alleles),
            "n_protein_annotations": sum(len(row.get("protein_annotations", []))
                                         for row in alpha_alleles),
            "jaspar_status": status_value(jaspar),
            "n_jaspar_alleles": len(jaspar_alleles),
            "n_motif_hits": sum(len(row.get("motif_hits", [])) for row in jaspar_alleles),
            "clinvar_status": status_value(clinvar),
            "n_clinvar_links": count_nested(clinvar, "links"),
            "n_clinvar_direct_records": count_nested(clinvar, "direct_records"),
            "n_clinvar_related_set_records": count_nested(clinvar, "related_set_records"),
            "gwas_catalog_status": status_value(gwas),
            "n_gwas_associations": count_nested(gwas, "associations"),
            "gnomad_status": status_value(gnomad),
            "n_gnomad_variant_links": len(gnomad_links),
            "n_gnomad_verified_variants": sum(item.get("variant") is not None
                                               for item in gnomad_links),
            "gtr_status": status_value(gtr),
            "n_gtr_test_matches": count_nested(gtr, "test_matches"),
            "alpha_fold_json": alpha, "jaspar_json": jaspar,
            "clinvar_json": clinvar, "gwas_catalog_json": gwas,
            "gnomad_json": gnomad, "gtr_json": gtr,
        })
    return rows


def write_sqlite(path: Path, records: list[dict], manifest: dict) -> None:
    partial = path.with_suffix(path.suffix + ".part")
    partial.unlink(missing_ok=True)
    database = sqlite3.connect(partial)
    try:
        database.execute("PRAGMA foreign_keys=ON")
        database.execute(
            "CREATE TABLE combined_data ("
            "rank INTEGER NOT NULL UNIQUE, rsid TEXT PRIMARY KEY, "
            "current_rsid TEXT, was_merged INTEGER NOT NULL, "
            "alpha_fold_json TEXT, jaspar_json TEXT, clinvar_json TEXT, "
            "gwas_catalog_json TEXT, gnomad_json TEXT, gtr_json TEXT, "
            "record_json TEXT NOT NULL)"
        )
        database.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT NOT NULL)")
        fields = list(SOURCE_NAMES)
        for record in records:
            dumped = [json.dumps(record.get(field), ensure_ascii=False, sort_keys=True,
                                 separators=(",", ":")) if record.get(field) is not None else None
                      for field in fields]
            database.execute(
                "INSERT INTO combined_data VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (record["rank"], record["rsid"], record["identity"]["current_rsid"],
                 int(record["identity"]["was_merged"]), *dumped,
                 json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"))),
            )
        database.execute(
            "INSERT INTO metadata VALUES (?,?)",
            ("run_manifest", json.dumps(manifest, ensure_ascii=False, sort_keys=True)),
        )
        integrity = database.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            raise RuntimeError(f"SQLite integrity check failed: {integrity}")
        database.commit()
    except Exception:
        database.close()
        partial.unlink(missing_ok=True)
        raise
    else:
        database.close()
        partial.replace(path)


def write_readme(path: Path, manifest: dict) -> None:
    sources = ", ".join(manifest["present_sources"])
    missing = ", ".join(manifest["missing_sources"]) or "none"
    text = f"""# Combined six-source rsID dataset

This folder combines the normalized AlphaFold, JASPAR, ClinVar, GWAS Catalog,
gnomAD, and GTR generator outputs without making cross-source biological
claims. There is exactly one top-level record per supplied rsID.

- Input rsIDs: {manifest['n_input_rsids']}
- Present sources: {sources}
- Missing sources: {missing}

`combined_data.json` is the authoritative nested representation.
`combined_summary.tsv` has one row per rsID and keeps every source object in a
JSON-valued column, alongside useful counts. `combined_data.sqlite` mirrors
those rows with one JSON column per source. `run_manifest.json` records all
consumed files, SHA-256 hashes, embedded source manifests, and source-reported
errors.

One-to-many records are nested by their actual relation keys: genomic variant,
ClinVar record, GWAS source-row digest, gnomAD variant/callset, or GTR test.
They are never flattened into a Cartesian product. Source-specific cautions
and evidence scopes remain authoritative.
"""
    (path / "README.md").write_text(text, encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    rsids = select_panel_rsids(args.input, require_exact=False)
    dataset_root = Path(args.dataset_root)
    out = Path(args.out) if args.out else dataset_root / "combined_dataset"
    out.mkdir(parents=True, exist_ok=True)
    # A failed rerun must not leave prior combined artifacts looking current.
    for filename in OUTPUT_FILENAMES:
        (out / filename).unlink(missing_ok=True)
    input_path = Path(args.input)
    source_records: dict[str, dict[str, dict]] = {}
    provenance: dict[str, dict] = {}
    missing_sources = []
    for source in SOURCE_NAMES:
        directory = dataset_root / SOURCE_DIRS[source]
        if not directory.is_dir():
            missing_sources.append(source)
            continue
        records, source_meta = LOADERS[source](directory, dataset_root, rsids)
        source_records[source] = records
        provenance[source] = source_meta
    if missing_sources and not args.allow_partial:
        raise FileNotFoundError(
            "All six dataset directories are required; missing "
            f"{missing_sources}. Use --allow-partial to combine available sources."
        )
    if not source_records:
        raise FileNotFoundError("No dataset source directories were found to combine")

    combined = []
    for rank, rsid in enumerate(rsids, 1):
        by_source = {source: source_records.get(source, {}).get(rsid)
                     for source in SOURCE_NAMES}
        identity_values: dict[str, str] = {}
        merged_values: dict[str, bool] = {}
        for source in ("clinvar", "gwas_catalog", "gnomad", "gtr"):
            value = by_source.get(source)
            if value:
                resolution = value.get("dbsnp_resolution", {})
                current = str(resolution.get("current_rsid", ""))
                if current:
                    identity_values[source] = current
                merged_values[source] = bool(resolution.get("was_merged", False))
        unique_current = sorted(set(identity_values.values()))
        if len(unique_current) > 1:
            raise ValueError(
                f"Conflicting current dbSNP identifiers for {rsid}: {identity_values}"
            )
        combined.append({
            "rank": rank,
            "rsid": rsid,
            "identity": {
                "current_rsid": unique_current[0] if unique_current else rsid,
                "was_merged": any(merged_values.values()),
                "current_rsid_by_source": identity_values,
                "was_merged_by_source": merged_values,
            },
            **by_source,
        })

    summary = build_summary(combined)
    preliminary_manifest = {
        "schema_version": "1.0",
        "generator": Path(__file__).name,
        "finished_at_utc": utc_now(),
        "input": str(input_path.resolve()),
        "input_sha256": sha256_file(input_path),
        "dataset_root": str(dataset_root.resolve()),
        "output": str(out.resolve()),
        "allow_partial": bool(args.allow_partial),
        "n_input_rsids": len(rsids),
        "input_rsids": rsids,
        "present_sources": [source for source in SOURCE_NAMES if source in source_records],
        "missing_sources": missing_sources,
        "source_provenance": provenance,
        "output_files": {},
    }
    write_json(out / "combined_data.json", combined)
    write_tsv(out / "combined_summary.tsv", summary, SUMMARY_FIELDS)
    write_sqlite(out / "combined_data.sqlite", combined, preliminary_manifest)
    preliminary_manifest["output_files"] = {
        filename: {
            "size_bytes": (out / filename).stat().st_size,
            "sha256": sha256_file(out / filename),
        }
        for filename in ("combined_data.json", "combined_summary.tsv", "combined_data.sqlite")
    }
    write_json(out / "run_manifest.json", preliminary_manifest)
    write_readme(out, preliminary_manifest)
    print(f"Combined {len(rsids)} rsIDs from {len(source_records)} sources. See {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT),
                        help="rsID TSV or one-rsID-per-line input")
    parser.add_argument("--dataset-root", default=str(DEFAULT_DATASET_ROOT),
                        help="Parent containing the six *_dataset directories")
    parser.add_argument("--out", help="Output directory (default: <dataset-root>/combined_dataset)")
    parser.add_argument("--allow-partial", action="store_true",
                        help="Allow missing source directories; present sources must still be complete")
    return parser


def main() -> int:
    try:
        return run(build_parser().parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
