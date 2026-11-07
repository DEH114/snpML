#!/usr/bin/env Rscript
# =============================================================================
# snp_evidence_ranker.R
#
# Rank the quality of evidence for a single SNP (dbSNP rsID) using six sources:
#   ClinVar, gnomAD, GWAS Catalog, AlphaFold DB (+ AlphaMissense), JASPAR, NIH GTR
#
# Output (JSON + console report):
#   confidence        0-100 score with band: High (>=80), Moderate (50-79),
#                     Low (25-49), Emerging (<25)
#   evidence_summary  items tagged direct / strong / indirect / mixed / emerging
#   limitations       enumerated per query
#   audit             one entry per reference, flagged for citation integrity
#   conflicts         each detected disagreement, classified as resolved by a
#                     newer reviewed correction (booster) or unresolved (penalty)
#
# Conflict model: concordant evidence is scored at baseline. When a newer,
# equal-or-higher-tier source with a review signal contradicts and supersedes
# an older one, the superseded assertions lose their weight and a capped
# "correction booster" is added on top of the baseline score.
#
# Usage:
#   Rscript snp_evidence_ranker.R rs429358
#   Rscript snp_evidence_ranker.R rs429358 --out=rs429358.json
#   Rscript snp_evidence_ranker.R rs113488022 --sources=clinvar,gnomad,alphafold --quiet
#
# From an R session:
#   source("snp_evidence_ranker.R")
#   res <- rank_snp_evidence("rs429358", out = "rs429358.json")
#
# Requires: httr2 (>= 1.0.0), jsonlite
# Optional: xml2 (enables submission-level ClinVar conflict analysis)
# Optional environment variables (raise NCBI E-utilities rate limits):
#   NCBI_API_KEY, NCBI_EMAIL
# =============================================================================

for (.pkg in c("httr2", "jsonlite")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with install.packages('%s').",
                 .pkg, .pkg), call. = FALSE)
  }
}
rm(.pkg)

TOOL_VERSION <- "1.1.0"

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
CONFIG <- list(
  user_agent      = paste0("snp-evidence-ranker/", TOOL_VERSION, " (R; httr2)"),
  timeout_s       = 60,
  ncbi_api_key    = Sys.getenv("NCBI_API_KEY", ""),
  ncbi_email      = Sys.getenv("NCBI_EMAIL", ""),

  ensembl_api     = "https://rest.ensembl.org",
  eutils_api      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
  gnomad_api      = "https://gnomad.broadinstitute.org/api",
  gnomad_dataset  = "gnomad_r4",
  gwas_v1_api     = "https://www.ebi.ac.uk/gwas/rest/api",
  gwas_v2_api     = "https://www.ebi.ac.uk/gwas/rest/api/v2",
  alphafold_api   = "https://alphafold.ebi.ac.uk/api",
  jaspar_release  = "JASPAR2024",
  jaspar_pfm_url  = paste0("https://jaspar.elixir.no/download/data/2024/CORE/",
                           "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"),
  cache_dir       = tools::R_user_dir("snp_evidence_ranker", which = "cache"),

  # Retrieval limits
  clinvar_max_records    = 15,
  clinvar_max_pmids      = 20,
  gwas_sig_p             = 5e-8,
  gwas_max_trait_items   = 12,
  gwas_max_study_lookups = 60,
  gtr_max_records        = 100,

  # Interpretation thresholds
  gnomad_min_an       = 10000,  # total allele number for "robust" frequency
  gnomad_pop_min_an   = 2000,   # min AN for a population to count toward max AF
  gwas_small_n        = 2000,
  gwas_large_n        = 50000,
  jaspar_flank        = 30,
  jaspar_min_rel      = 0.80,   # relative PWM score to call a motif match
  jaspar_min_delta    = 0.10,   # min |alt - ref| relative score change
  jaspar_max_items    = 5,
  stale_years         = 5,      # ClinVar classification older than this is flagged

  # Conflict resolution (review-and-correction model)
  gnomad_release_date     = "2023-11-01",  # release date of the gnomAD dataset above
  correction_min_gap_days = 180,  # corrector must postdate every opposing assertion by this much
  gwas_correction_n_ratio = 2,    # later GWAS must be >= this many times larger (implicit review)
  booster_explicit        = 10,   # score points: expert panel / guideline / ClinVar flag / meta-analysis
  booster_implicit        = 5,    # score points: >= 2 independent newer sources, or a much larger study
  booster_decay           = 0.5,  # i-th correction counts booster * decay^(i-1)
  booster_cap             = 15,   # maximum total booster

  # Scoring model
  strength_points     = c(direct = 30, strong = 18, indirect = 8, mixed = 6, emerging = 3),
  within_source_decay = 0.6,    # diminishing returns for repeated items from one source
  saturation_k        = 45      # score = 100 * (1 - exp(-points / k))
)

SOURCES <- c(clinvar = "ClinVar", gnomad = "gnomAD", gwas = "GWAS Catalog",
             alphafold = "AlphaFold", jaspar = "JASPAR", gtr = "NIH GTR")

STRENGTH_DEFS <- list(
  direct   = "Primary studies or curated assertions directly answering the query",
  strong   = "Robust indirect or mechanistic evidence",
  indirect = "Related but not on-point evidence",
  mixed    = "Conflicting results across studies or sources",
  emerging = "Early-stage or hypothesis-generating data"
)

BAND_DEFS <- list(
  High     = "Score >= 80: multiple direct, replicated studies (RCT/meta-analysis-grade where available)",
  Moderate = "Score 50-79: good but heterogeneous or limited evidence",
  Low      = "Score 25-49: sparse, indirect, or preclinical/computational-only data",
  Emerging = "Score < 25: very early or conflicting signals"
)

ACC_PATTERNS <- list(
  "ClinVar"      = "^VCV[0-9]{9}(\\.[0-9]+)?$|^uid[0-9]+$",
  "GWAS Catalog" = "^GCST[0-9]+$",
  "AlphaFold DB" = "^AF-[A-Z0-9]+-F[0-9]+$",
  "JASPAR"       = "^MA[0-9]{4}\\.[0-9]+$",
  "NIH GTR"      = "^GTR[0-9]+(\\.[0-9]+)?$|^uid[0-9]+$",
  "gnomAD"       = "^(chr)?[0-9XYM]+-[0-9]+-[ACGT]+-[ACGT]+$|^rs[0-9]+$"
)

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

pluck <- function(x, ...) {
  for (k in list(...)) {
    if (is.null(x)) return(NULL)
    x <- tryCatch(x[[k]], error = function(e) NULL)
  }
  x
}

get_any <- function(x, keys) {
  if (!is.list(x)) return(NULL)
  for (k in keys) {
    v <- x[[k]]
    if (!is.null(v) && length(v) > 0) return(v)
  }
  NULL
}

chr1 <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  v <- unlist(x, use.names = FALSE)
  if (!length(v) || is.na(v[1])) return(default)
  as.character(v[1])
}

num1 <- function(x, default = NA_real_) {
  v <- suppressWarnings(as.numeric(chr1(x, NA_character_)))
  if (is.na(v)) default else v
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else max(x)
}

fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format = "e", digits = 1)

year_of <- function(x) {
  y <- suppressWarnings(as.integer(substr(chr1(x, ""), 1, 4)))
  if (is.na(y)) NA_integer_ else y
}

# -----------------------------------------------------------------------------
# HTTP layer (retries, per-host throttling, never throws)
# -----------------------------------------------------------------------------
http_request <- function(url, query = NULL, body = NULL, realm = NULL, rate = 5,
                         parse = c("json", "text")) {
  parse <- match.arg(parse)
  req <- httr2::request(url)
  req <- httr2::req_user_agent(req, CONFIG$user_agent)
  req <- httr2::req_timeout(req, CONFIG$timeout_s)
  req <- httr2::req_retry(req, max_tries = 4,
                          is_transient = function(resp)
                            httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504))
  req <- httr2::req_error(req, is_error = function(resp) FALSE)
  if (length(query)) req <- do.call(httr2::req_url_query, c(list(req), query))
  if (!is.null(realm)) req <- httr2::req_throttle(req, rate = rate, realm = realm)
  if (parse == "json") req <- httr2::req_headers(req, Accept = "application/json")
  if (!is.null(body)) req <- httr2::req_body_json(req, body)

  resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
  if (inherits(resp, "error")) {
    return(list(ok = FALSE, status = NA_integer_, data = NULL,
                error = conditionMessage(resp), url = url))
  }
  status <- httr2::resp_status(resp)
  txt <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
  data <- if (parse == "json") {
    tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
  } else txt
  ok <- status >= 200 && status < 300 && !is.null(data)
  list(ok = ok, status = status, data = data,
       error = if (!ok) sprintf("HTTP %s from %s", status, url) else NULL, url = url)
}

ensembl_get <- function(path, query = list()) {
  http_request(paste0(CONFIG$ensembl_api, path),
               query = c(list(`content-type` = "application/json"), query),
               realm = "ensembl", rate = 10)
}

eutils <- function(util, query) {
  query$retmode <- "json"
  query$tool <- "snp_evidence_ranker"
  if (nzchar(CONFIG$ncbi_api_key)) query$api_key <- CONFIG$ncbi_api_key
  if (nzchar(CONFIG$ncbi_email)) query$email <- CONFIG$ncbi_email
  rate <- if (nzchar(CONFIG$ncbi_api_key)) 9 else 2.5
  http_request(sprintf("%s/%s.fcgi", CONFIG$eutils_api, util),
               query = query, realm = "ncbi", rate = rate)
}

# -----------------------------------------------------------------------------
# Core data structures
# -----------------------------------------------------------------------------
CONFLICT_CONTEXTS <- c("baseline_no_conflict", "corrected_by_review", "unresolved_conflict")

new_item <- function(source, claim, strength, quality = 1, replicated = FALSE,
                     study_ids = character(), direction = 0L, details = list(),
                     refs = character(), conflict_context = "baseline_no_conflict",
                     conflict_key = NA_character_) {
  strength <- match.arg(strength, names(CONFIG$strength_points))
  conflict_context <- match.arg(conflict_context, CONFLICT_CONTEXTS)
  list(source = source, claim = claim, strength = strength,
       quality = max(0, min(1.5, quality)), replicated = isTRUE(replicated),
       study_ids = unique(as.character(study_ids)), direction = as.integer(direction),
       details = details, refs = unique(as.character(refs)),
       conflict_context = conflict_context, conflict_key = conflict_key)
}

new_ref <- function(raw_id, ref_type = c("pmid", "database_record"), source,
                    database = NULL, claimed = list(), flags = character(),
                    url = NA_character_) {
  ref_type <- match.arg(ref_type)
  raw_id <- as.character(raw_id)
  db <- if (!is.null(database)) database else if (ref_type == "pmid") "PubMed" else source
  ref_id <- if (ref_type == "pmid") paste0("PMID:", raw_id) else paste0(db, ":", raw_id)
  if (ref_type == "pmid" && is.na(url)) url <- sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", raw_id)
  list(ref_id = ref_id, raw_id = raw_id, ref_type = ref_type, database = db,
       cited_by = source, claimed = claimed, flags = flags, url = url)
}

module_result <- function(source, status = c("ok", "no_data", "skipped", "error"),
                          items = list(), refs = list(), notes = character(),
                          flags = list(), raw = list(), error = NULL) {
  list(source = source, status = match.arg(status), items = items, refs = refs,
       notes = notes, flags = flags, raw = raw, error = error)
}

# -----------------------------------------------------------------------------
# Conflict resolution: review-and-correction model
#
# An "assertion" is one dated position on one question (a ClinVar submission's
# classification of a VCV record, or one GWAS study's effect direction for a
# trait). A disagreement counts as RESOLVED only when one position:
#   (a) postdates every opposing live assertion by >= correction_min_gap_days,
#   (b) comes from an equal-or-higher review tier than the opposing assertions,
#   (c) carries a review signal:
#         explicit - expert panel / practice guideline, ClinVar-flagged opposing
#                    submissions, or a later meta-analysis
#         implicit - >= 2 independent newer submitters/studies, or a later study
#                    >= gwas_correction_n_ratio times larger
# Newer alone is never enough. Anything else is UNRESOLVED.
# -----------------------------------------------------------------------------
safe_date <- function(x) {
  if (inherits(x, "Date")) return(x[1])
  s <- substr(chr1(x, ""), 1, 10)
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", s)) return(as.Date(NA))
  tryCatch(as.Date(s), error = function(e) as.Date(NA))
}

new_assertion <- function(id, source, group, category, date, tier,
                          submitter = NA_character_, flagged = FALSE, n = NA_real_,
                          meta = FALSE, label = NA_character_) {
  list(id = as.character(id), source = source, group = group, category = category,
       date = safe_date(date), tier = as.numeric(tier), submitter = submitter,
       flagged = isTRUE(flagged), n = n, meta = isTRUE(meta), label = label)
}

resolve_conflict <- function(asr, source, group) {
  asr <- Filter(function(a) !is.na(a$category), asr)
  summ <- function(a) list(id = a$id, position = a$category,
                           date = if (is.na(a$date)) NA_character_ else format(a$date),
                           review_tier = a$tier, submitter = a$submitter, flagged = a$flagged,
                           sample_size = a$n, label = a$label)
  cats <- unique(vapply(asr, `[[`, "", "category"))
  base <- list(source = source, group = group, positions = as.list(cats))
  if (length(cats) < 2) return(c(base, list(status = "concordant")))
  tier_of <- function(x) vapply(x, `[[`, numeric(1), "tier")

  for (cc in cats) {
    new_side <- Filter(function(a) a$category == cc && !a$flagged, asr)
    old_side <- Filter(function(a) a$category != cc, asr)
    if (!length(new_side)) next
    old_live <- Filter(function(a) !a$flagged, old_side)
    any_flagged_old <- any(vapply(old_side, `[[`, logical(1), "flagged"))
    signal <- NA_character_; basis <- character()

    if (!length(old_live)) {
      correctors <- new_side
      signal <- "explicit"
      basis <- "every opposing submission has been flagged by ClinVar as outdated or erroneous"
    } else {
      old_dates <- do.call(c, lapply(old_live, `[[`, "date"))
      if (anyNA(old_dates)) next                      # cannot order undated claims
      latest_old <- max(old_dates)
      correctors <- Filter(function(a) !is.na(a$date) &&
                             as.numeric(a$date - latest_old) >= CONFIG$correction_min_gap_days,
                           new_side)
      if (!length(correctors)) next                   # (a) not newer
      if (max(tier_of(correctors)) < max(tier_of(old_live))) next   # (b) lower authority
      subs <- unique(stats::na.omit(vapply(correctors, function(a) chr1(a$submitter), "")))
      old_n <- safe_max(vapply(old_live, function(a) num1(a$n), numeric(1)))
      new_n <- safe_max(vapply(correctors, function(a) num1(a$n), numeric(1)))
      if (source == "ClinVar" && max(tier_of(correctors)) >= 3) {
        signal <- "explicit"
        basis <- c(basis, "an expert-panel or practice-guideline review postdates the opposing submissions")
      }
      if (any(vapply(correctors, `[[`, logical(1), "meta"))) {
        signal <- "explicit"
        basis <- c(basis, "a later meta-analysis re-examined the association")
      }
      if (any_flagged_old) {
        signal <- "explicit"
        basis <- c(basis, "some opposing submissions were flagged by ClinVar")
      }
      if (is.na(signal) && length(subs) >= 2) {
        signal <- "implicit"
        basis <- c(basis, sprintf("%d independent newer submitters/studies converged on this position", length(subs)))
      }
      if (is.na(signal) && !is.na(old_n) && !is.na(new_n) &&
          new_n >= CONFIG$gwas_correction_n_ratio * old_n) {
        signal <- "implicit"
        basis <- c(basis, sprintf("a later study was %.1fx larger (N %s vs %s)", new_n / old_n,
                                  format(new_n, big.mark = ","), format(old_n, big.mark = ",")))
      }
      if (is.na(signal)) next                         # (c) newer but no evidence of review
    }

    subs_all <- unique(stats::na.omit(vapply(correctors, function(a) chr1(a$submitter), "")))
    cdates <- do.call(c, lapply(correctors, `[[`, "date"))
    return(c(base, list(
      status = "resolved_correction", current_position = cc,
      superseded_positions = as.list(setdiff(cats, cc)),
      review_signal = signal, signal_basis = as.list(basis),
      max_corrector_tier = max(tier_of(correctors)),
      n_corrector_sources = length(subs_all),
      corrector_earliest_date = if (length(cdates) && any(!is.na(cdates)))
        format(min(cdates, na.rm = TRUE)) else NA_character_,
      correctors = lapply(correctors, summ), superseded = lapply(old_side, summ))))
  }
  c(base, list(status = "unresolved",
               reason = "No position postdates all opposing assertions with an equal-or-higher review tier and a review signal.",
               assertions = lapply(asr, summ)))
}

# -----------------------------------------------------------------------------
# 0. Resolve the SNP (Ensembl variation + VEP)
# -----------------------------------------------------------------------------
resolve_variant <- function(rsid) {
  r <- ensembl_get(paste0("/variation/human/", rsid))
  if (!r$ok) {
    stop(sprintf("Could not resolve %s via Ensembl REST (%s).", rsid,
                 r$error %||% "unknown error"), call. = FALSE)
  }
  d <- r$data
  std_chr <- c(as.character(1:22), "X", "Y", "MT")
  maps <- Filter(function(m) identical(chr1(m$assembly_name), "GRCh38") &&
                   chr1(m$seq_region_name, "") %in% std_chr,
                 d$mappings %||% list())
  if (!length(maps)) stop(sprintf("%s has no GRCh38 primary-assembly mapping.", rsid), call. = FALSE)
  m <- maps[[1]]

  alleles <- strsplit(chr1(m$allele_string, ""), "/", fixed = TRUE)[[1]]
  if (num1(m$strand, 1) == -1) {
    alleles <- vapply(alleles, function(a) {
      if (a == "-") a else paste(rev(strsplit(chartr("ACGT", "TGCA", a), "")[[1]]), collapse = "")
    }, "", USE.NAMES = FALSE)
  }
  ref <- alleles[1]
  alts <- alleles[-1]

  var <- list(
    rsid = rsid, rsid_current = chr1(d$name, rsid),
    synonyms = as.character(unlist(d$synonyms)),
    var_class = chr1(d$var_class), most_severe_consequence = chr1(d$most_severe_consequence),
    assembly = "GRCh38", chr = chr1(m$seq_region_name),
    start = num1(m$start), end = num1(m$end),
    ref = ref, alts = alts, allele_string = paste(alleles, collapse = "/"),
    is_snv = nchar(ref) == 1 && length(alts) > 0 && all(nchar(alts) == 1) &&
      ref != "-" && all(alts != "-"),
    n_grch38_mappings = length(maps),
    minor_allele = chr1(d$minor_allele), global_maf = num1(d$MAF),
    gene = NA_character_, gene_id = NA_character_, transcript_id = NA_character_,
    consequence_terms = character(), protein_pos = NA_real_, ref_aa = NA_character_,
    uniprot = NA_character_, missense = list(),
    regulatory_overlaps = 0L, motif_overlaps = 0L, vep_ok = FALSE
  )

  v <- ensembl_get(paste0("/vep/human/id/", var$rsid_current),
                   list(canonical = 1, uniprot = 1, mane = 1))
  if (v$ok && length(v$data)) {
    var$vep_ok <- TRUE
    vd <- v$data[[1]]
    var$regulatory_overlaps <- length(vd$regulatory_feature_consequences %||% list())
    var$motif_overlaps <- length(vd$motif_feature_consequences %||% list())
    tcs <- vd$transcript_consequences %||% list()
    if (length(tcs)) {
      is_canon <- vapply(tcs, function(t) num1(t$canonical, 0) == 1, logical(1))
      is_mane  <- vapply(tcs, function(t) !is.null(t$mane_select), logical(1))
      has_prot <- vapply(tcs, function(t) !is.null(t$protein_start), logical(1))
      pick <- tcs[[which.max(4 * is_mane + 2 * is_canon + has_prot)]]

      var$gene <- chr1(pick$gene_symbol)
      var$gene_id <- chr1(pick$gene_id)
      var$transcript_id <- chr1(pick$transcript_id)
      var$consequence_terms <- as.character(unlist(pick$consequence_terms))
      var$protein_pos <- num1(pick$protein_start)
      if (!is.null(pick$amino_acids)) {
        var$ref_aa <- substr(strsplit(chr1(pick$amino_acids), "/", fixed = TRUE)[[1]][1], 1, 1)
      }
      sp <- as.character(unlist(pick$swissprot))
      if (length(sp)) var$uniprot <- sub("[-.].*$", "", sp[1])

      same_tx <- Filter(function(t) identical(chr1(t$transcript_id), var$transcript_id), tcs)
      for (t in same_tx) {
        if ("missense_variant" %in% unlist(t$consequence_terms) && !is.null(t$amino_acids)) {
          aa <- strsplit(chr1(t$amino_acids), "/", fixed = TRUE)[[1]]
          pos <- num1(t$protein_start)
          var$missense[[length(var$missense) + 1]] <- list(
            allele = chr1(t$variant_allele), pos = pos, ref_aa = aa[1], alt_aa = aa[2],
            protein_change = paste0(aa[1], pos, aa[2]))
        }
      }
    }
  }
  var
}

# -----------------------------------------------------------------------------
# 1. ClinVar (NCBI E-utilities)
# -----------------------------------------------------------------------------
review_stars <- function(review) {
  r <- tolower(review %||% "")
  if (grepl("practice guideline", r)) return(4L)
  if (grepl("expert panel", r)) return(3L)
  if (grepl("multiple submitters", r) && grepl("no conflict", r)) return(2L)
  if (grepl("conflicting", r)) return(1L)
  if (grepl("single submitter", r)) return(1L)
  0L
}

clin_direction <- function(desc) {
  d <- tolower(desc %||% "")
  if (grepl("conflicting", d)) return(0L)
  if (grepl("pathogenic", d) && !grepl("benign", d)) return(1L)
  if (grepl("benign", d) && !grepl("pathogenic", d)) return(-1L)
  0L
}

clinvar_pmids <- function(uid) {
  el <- eutils("elink", list(dbfrom = "clinvar", db = "pubmed", id = uid,
                             linkname = "clinvar_pubmed"))
  if (!el$ok) return(character())
  dbs <- pluck(el$data, "linksets", 1, "linksetdbs") %||% list()
  head(unique(as.character(unlist(lapply(dbs, function(d) d$links)))), CONFIG$clinvar_max_pmids)
}

clin_category <- function(desc) {
  d <- tolower(desc %||% "")
  if (!nzchar(d) || grepl("conflicting", d)) return(NA_character_)
  if (grepl("pathogenic", d) && !grepl("benign", d)) return("pathogenic")
  if (grepl("benign", d) && !grepl("pathogenic", d)) return("benign")
  if (grepl("uncertain", d)) return("uncertain")
  NA_character_   # risk factor, drug response, other: outside the conflict model
}

cat_direction <- function(cat) {
  if (is.na(cat)) return(0L)
  switch(cat, pathogenic = 1L, benign = -1L, 0L)
}

scv_tier <- function(review) {
  r <- tolower(review %||% "")
  if (grepl("practice guideline", r)) 4 else if (grepl("expert panel", r)) 3
  else if (grepl("criteria provided", r)) 1 else 0
}

# Submission-level (SCV) records for one ClinVar variation, as dated assertions.
# Returns NULL when xml2 is missing or the fetch fails.
clinvar_scvs <- function(uid, vcv) {
  if (!requireNamespace("xml2", quietly = TRUE)) return(NULL)
  url <- sprintf("%s/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=%s&tool=snp_evidence_ranker%s%s",
                 CONFIG$eutils_api, uid,
                 if (nzchar(CONFIG$ncbi_api_key)) paste0("&api_key=", CONFIG$ncbi_api_key) else "",
                 if (nzchar(CONFIG$ncbi_email)) paste0("&email=", utils::URLencode(CONFIG$ncbi_email, reserved = TRUE)) else "")
  r <- http_request(url, parse = "text", realm = "ncbi",
                    rate = if (nzchar(CONFIG$ncbi_api_key)) 9 else 2.5)
  if (!r$ok) return(NULL)
  doc <- tryCatch(xml2::read_xml(r$data), error = function(e) NULL)
  if (is.null(doc)) return(NULL)
  cas <- xml2::xml_find_all(doc, ".//ClinicalAssertion")
  out <- list()
  for (k in seq_along(cas)) {
    ca <- cas[[k]]
    acc <- xml2::xml_find_first(ca, "./ClinVarAccession")
    cls_node <- xml2::xml_find_first(ca, "./Classification | ./Interpretation")
    desc <- xml2::xml_text(xml2::xml_find_first(cls_node, "./GermlineClassification | ./Description"))
    category <- clin_category(desc)
    if (is.na(category)) next                       # somatic / risk factor / other
    review <- xml2::xml_text(xml2::xml_find_first(
      ca, "./ReviewStatus | ./Classification/ReviewStatus | ./Interpretation/ReviewStatus"))
    flagged <- grepl("flagged", tolower(chr1(review, ""))) ||
      length(xml2::xml_find_all(ca, ".//*[contains(local-name(), 'Flagged')]")) > 0
    date <- safe_date(xml2::xml_attr(cls_node, "DateLastEvaluated"))
    if (is.na(date)) date <- safe_date(xml2::xml_attr(ca, "SubmissionDate"))
    if (is.na(date)) date <- safe_date(xml2::xml_attr(acc, "DateUpdated"))
    submitter <- xml2::xml_attr(acc, "SubmitterName")
    out[[length(out) + 1]] <- new_assertion(
      id = chr1(xml2::xml_attr(acc, "Accession"), sprintf("%s-scv%d", vcv, k)),
      source = "ClinVar", group = vcv, category = category, date = date,
      tier = if (flagged) 0 else scv_tier(review),
      submitter = chr1(xml2::xml_attr(acc, "OrgID"), chr1(submitter)),
      flagged = flagged,
      label = sprintf("%s: %s (%s)", chr1(submitter, "unknown submitter"), desc,
                      chr1(review, "no review status")))
  }
  out
}

query_clinvar <- function(var) {
  src <- "ClinVar"
  s <- eutils("esearch", list(db = "clinvar", term = var$rsid_current, retmax = 100))
  if (!s$ok) return(module_result(src, "error", error = s$error))
  ids <- as.character(unlist(pluck(s$data, "esearchresult", "idlist")))
  if (!length(ids)) {
    return(module_result(src, "no_data", notes = sprintf("No ClinVar records for %s.", var$rsid_current)))
  }
  truncated <- length(ids) > CONFIG$clinvar_max_records
  ids <- head(ids, CONFIG$clinvar_max_records)

  sm <- eutils("esummary", list(db = "clinvar", id = paste(ids, collapse = ",")))
  if (!sm$ok) return(module_result(src, "error", error = sm$error))
  res <- pluck(sm$data, "result")
  rs_num <- sub("^rs", "", var$rsid_current)
  this_year <- as.integer(format(Sys.Date(), "%Y"))

  items <- list(); refs <- list(); recs <- list(); notes <- character(); conflicts <- list()
  single_submitter <- FALSE; stale <- character(); unconfirmed <- 0L; scv_failed <- 0L

  for (uid in ids) {
    rec <- res[[uid]]
    if (is.null(rec) || !is.null(rec$error)) next

    xrefs <- unlist(lapply(rec$variation_set %||% list(), function(vs)
      vapply(vs$variation_xrefs %||% list(), function(x)
        paste0(chr1(x$db_source, ""), ":", chr1(x$db_id, "")), "")))
    xref_ok <- if (length(xrefs)) any(xrefs == paste0("dbSNP:", rs_num)) else NA
    if (isFALSE(xref_ok)) next   # text match, but a different variant (e.g. haplotype)

    gc <- rec$germline_classification
    if (is.null(gc) || !nzchar(chr1(gc$description, ""))) gc <- rec$clinical_significance
    desc <- chr1(gc$description, "")
    if (!nzchar(desc)) {
      notes <- c(notes, sprintf("ClinVar uid %s has no germline classification (somatic/oncogenicity only); not scored.", uid))
      next
    }
    acc <- chr1(rec$accession, paste0("uid", uid))
    review <- chr1(gc$review_status, "")
    last_eval <- chr1(gc$last_evaluated)
    traits <- unique(stats::na.omit(vapply(gc$trait_set %||% rec$trait_set %||% list(),
                                            function(t) chr1(t$trait_name), "")))
    stars <- review_stars(review)
    dirn <- clin_direction(desc)
    conflict <- grepl("conflicting", tolower(paste(desc, review)))
    vus <- grepl("uncertain significance", tolower(desc)) && !conflict
    ly <- year_of(last_eval)
    if (!is.na(ly) && (this_year - ly) > CONFIG$stale_years) stale <- c(stale, acc)
    if (is.na(xref_ok)) unconfirmed <- unconfirmed + 1L

    pmids <- clinvar_pmids(uid)
    db_flags <- if (is.na(xref_ok)) "dbsnp_xref_unconfirmed" else character()
    refs <- c(refs, list(new_ref(acc, "database_record", src, "ClinVar", flags = db_flags,
                                 url = sprintf("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/", uid))))
    for (p in pmids) refs <- c(refs, list(new_ref(p, "pmid", src)))

    # Baseline strength from the aggregate (VCV-level) classification
    if (conflict)         { strength <- "mixed";    q <- 1.0; repl <- FALSE }
    else if (vus)         { strength <- "emerging"; q <- 0.8; repl <- FALSE }
    else if (stars >= 3)  { strength <- "direct";   q <- 1.5; repl <- TRUE  }
    else if (stars == 2)  { strength <- "direct";   q <- 1.2; repl <- TRUE  }
    else if (stars == 1)  { strength <- "direct";   q <- 0.7; repl <- FALSE; single_submitter <- TRUE }
    else                  { strength <- "emerging"; q <- 0.5; repl <- FALSE }

    trait_txt <- if (length(traits)) paste(head(traits, 3), collapse = "; ") else "unspecified condition"
    claim <- sprintf("%s classified '%s' for %s (review status: %s; %d star%s%s).",
                     acc, desc, trait_txt, if (nzchar(review)) review else "not provided",
                     stars, if (stars == 1) "" else "s",
                     if (!is.na(last_eval)) paste0("; last evaluated ", substr(last_eval, 1, 10)) else "")
    if (vus) claim <- paste(claim, "A VUS indicates insufficient evidence to classify.")

    # Submission-level conflict analysis
    scvs <- clinvar_scvs(uid, acc)
    if (is.null(scvs)) scv_failed <- scv_failed + 1L
    cf <- if (length(scvs)) resolve_conflict(scvs, src, acc) else NULL
    agg_cat <- clin_category(desc)
    ctx <- "baseline_no_conflict"

    if (!is.null(cf) && identical(cf$status, "resolved_correction")) {
      ctx <- "corrected_by_review"
      newcat <- cf$current_position
      sup_txt <- sprintf("%d older %s assertion(s) (%s)", length(cf$superseded),
                         paste(unlist(cf$superseded_positions), collapse = "/"),
                         paste(head(vapply(cf$superseded, `[[`, "", "id"), 5), collapse = ", "))
      basis_txt <- paste(unlist(cf$signal_basis), collapse = "; ")
      if (conflict || is.na(agg_cat) || agg_cat != newcat) {
        # The aggregate still reflects the superseded position: score the correction instead.
        strength <- if (newcat == "uncertain") "emerging" else "direct"
        q <- if (cf$max_corrector_tier >= 3) 1.5 else if (cf$n_corrector_sources >= 2) 1.2 else 0.9
        repl <- cf$max_corrector_tier >= 3 || cf$n_corrector_sources >= 2
        dirn <- cat_direction(newcat); conflict <- FALSE; vus <- newcat == "uncertain"
        claim <- sprintf("%s: ClinVar aggregate reads '%s', but %d submission(s) evaluated since %s classify the variant as %s for %s, superseding %s. Review signal: %s (%s).",
                         acc, desc, length(cf$correctors), cf$corrector_earliest_date, newcat,
                         trait_txt, sup_txt, cf$review_signal, basis_txt)
      } else {
        claim <- paste0(claim, sprintf(" This position superseded %s (review signal: %s; %s).",
                                       sup_txt, cf$review_signal, basis_txt))
      }
      conflicts[[length(conflicts) + 1]] <- cf
    } else if (!is.null(cf) && identical(cf$status, "unresolved")) {
      ctx <- "unresolved_conflict"
      conflicts[[length(conflicts) + 1]] <- cf
      if (conflict) {
        claim <- paste(claim, "No newer, equal-or-higher-tier reviewed submission supersedes the others.")
      } else {
        items[[length(items) + 1]] <- new_item(
          src, sprintf("Unresolved disagreement among ClinVar submissions for %s (positions: %s); no newer, equal-or-higher-tier reviewed submission supersedes the others.",
                       acc, paste(unlist(cf$positions), collapse = " vs ")),
          "mixed", 0.5, study_ids = paste0("clinvar-conflict:", acc),
          refs = paste0("ClinVar:", acc), conflict_context = "unresolved_conflict",
          conflict_key = paste(src, acc, sep = "|"))
      }
    }

    items[[length(items) + 1]] <- new_item(
      src, claim, strength, q, repl, study_ids = paste0("clinvar:", acc), direction = dirn,
      details = list(accession = acc, variation_id = uid, classification = desc,
                     review_status = review, stars = stars, last_evaluated = last_eval,
                     traits = as.list(traits), n_linked_pmids = length(pmids),
                     n_submissions_analysed = length(scvs %||% list()),
                     conflict_status = if (is.null(cf)) "not_assessed" else cf$status),
      refs = c(paste0("ClinVar:", acc), paste0("PMID:", pmids)),
      conflict_context = ctx,
      conflict_key = if (ctx == "baseline_no_conflict") NA_character_ else paste(src, acc, sep = "|"))

    plp_dates <- if (length(scvs)) do.call(c, lapply(Filter(function(a) a$category == "pathogenic" && !a$flagged, scvs), `[[`, "date")) else NULL
    latest_plp <- if (length(plp_dates) && any(!is.na(plp_dates))) format(max(plp_dates, na.rm = TRUE))
                  else if (dirn == 1L) substr(chr1(last_eval, ""), 1, 10) else NA_character_
    recs[[length(recs) + 1]] <- list(
      accession = acc,
      description = if (ctx == "corrected_by_review") paste(cf$current_position, "(after correction)") else desc,
      stars = stars, direction = dirn, conflict = conflict, vus = vus,
      latest_plp_date = latest_plp, conflict_status = if (is.null(cf)) "not_assessed" else cf$status)
  }

  if (!length(items)) {
    return(module_result(src, "no_data", notes = c(notes, "No germline ClinVar classifications matched this rsID.")))
  }
  if (scv_failed > 0) notes <- c(notes, sprintf("Submission-level records unavailable for %d ClinVar record(s)%s; their conflicts could not be assessed.",
                                                scv_failed, if (!requireNamespace("xml2", quietly = TRUE)) " (package xml2 not installed)" else ""))
  module_result(src, "ok", items, refs, notes,
                flags = list(truncated = truncated, single_submitter = single_submitter,
                             stale = stale, xref_unconfirmed = unconfirmed, scv_failed = scv_failed),
                raw = list(records = recs, conflicts = conflicts))
}

# -----------------------------------------------------------------------------
# 2. gnomAD (GraphQL)
# -----------------------------------------------------------------------------
GNOMAD_FIELDS <- "variant_id rsids
  exome  { ac an homozygote_count filters populations { id ac an } }
  genome { ac an homozygote_count filters populations { id ac an } }"

GNOMAD_Q_ID <- sprintf(
  "query Q($variantId: String!, $dataset: DatasetId!) { variant(variantId: $variantId, dataset: $dataset) { %s } }",
  GNOMAD_FIELDS)
GNOMAD_Q_RS <- sprintf(
  "query Q($rsid: String!, $dataset: DatasetId!) { variant(rsid: $rsid, dataset: $dataset) { %s } }",
  GNOMAD_FIELDS)

gnomad_fetch <- function(query, variables) {
  r <- http_request(CONFIG$gnomad_api, body = list(query = query, variables = variables),
                    realm = "gnomad", rate = 1)
  v <- pluck(r$data, "data", "variant")
  if (!is.null(v)) return(list(found = TRUE, v = v))
  errs <- pluck(r$data, "errors")
  msg <- if (length(errs)) paste(vapply(errs, function(e) chr1(e$message, ""), ""), collapse = "; ") else ""
  if (grepl("not found", msg, ignore.case = TRUE)) return(list(found = FALSE))
  if (isTRUE(r$ok) && !length(errs)) return(list(found = FALSE))
  list(found = NA, error = if (nzchar(msg)) msg else (r$error %||% "gnomAD request failed"))
}

summarise_gnomad <- function(v) {
  parts <- Filter(Negate(is.null), list(exome = v$exome, genome = v$genome))
  ac <- sum(vapply(parts, function(p) num1(p$ac, 0), numeric(1)))
  an <- sum(vapply(parts, function(p) num1(p$an, 0), numeric(1)))
  hom <- sum(vapply(parts, function(p) num1(p$homozygote_count, 0), numeric(1)))
  filters <- unique(as.character(unlist(lapply(parts, function(p) unlist(p$filters)))))
  pops <- list()
  for (p in parts) for (pp in p$populations %||% list()) {
    id <- chr1(pp$id, "")
    if (!nzchar(id) || grepl("_", id) || id %in% c("XX", "XY")) next
    cur <- pops[[id]] %||% c(ac = 0, an = 0)
    pops[[id]] <- cur + c(num1(pp$ac, 0), num1(pp$an, 0))
  }
  pop_af <- vapply(pops, function(x) if (x[2] >= CONFIG$gnomad_pop_min_an) x[1] / x[2] else NA_real_, numeric(1))
  max_pop <- if (length(pop_af) && any(!is.na(pop_af))) names(which.max(pop_af)) else NA_character_
  list(variant_id = chr1(v$variant_id), ac = ac, an = an, af = if (an > 0) ac / an else NA_real_,
       homozygotes = hom, filters = filters, max_pop = max_pop,
       max_pop_af = if (!is.na(max_pop)) unname(pop_af[[max_pop]]) else NA_real_,
       sequencing = names(parts))
}

query_gnomad <- function(var) {
  src <- "gnomAD"
  items <- list(); refs <- list(); notes <- character()
  low_an <- FALSE; filter_fail <- FALSE; max_af <- NA_real_; any_error <- NULL
  targets <- if (var$is_snv) {
    lapply(var$alts, function(a) list(alt = a, q = GNOMAD_Q_ID,
                                       vars = list(variantId = paste(var$chr, var$start, var$ref, a, sep = "-"),
                                                   dataset = CONFIG$gnomad_dataset)))
  } else {
    list(list(alt = paste(var$alts, collapse = ","), q = GNOMAD_Q_RS,
              vars = list(rsid = var$rsid_current, dataset = CONFIG$gnomad_dataset)))
  }

  for (tg in targets) {
    f <- gnomad_fetch(tg$q, tg$vars)
    vid <- tg$vars$variantId %||% tg$vars$rsid
    if (is.na(f$found)) { any_error <- f$error; next }
    if (!f$found) {
      items[[length(items) + 1]] <- new_item(
        src, sprintf("%s (alt %s) is absent from gnomAD %s. Absence is consistent with, but not proof of, a rare deleterious allele (ACMG PM2_Supporting-type evidence); it may also reflect low coverage.",
                     vid, tg$alt, CONFIG$gnomad_dataset),
        "indirect", 0.6, direction = 1L, study_ids = paste0("gnomad:", vid),
        details = list(variant_id = vid, dataset = CONFIG$gnomad_dataset, found = FALSE),
        refs = paste0("gnomAD:", vid))
      refs <- c(refs, list(new_ref(vid, "database_record", src, "gnomAD",
                                   url = sprintf("https://gnomad.broadinstitute.org/variant/%s?dataset=%s", vid, CONFIG$gnomad_dataset))))
      next
    }
    s <- summarise_gnomad(f$v)
    if (!is.na(s$variant_id)) vid <- s$variant_id
    af <- s$af; max_af <- safe_max(c(max_af, af, s$max_pop_af))
    pass <- !length(s$filters)
    if (!pass) filter_fail <- TRUE

    if (is.na(af) || s$ac == 0) {
      strength <- "indirect"; q <- 0.6; dirn <- 1L
      txt <- "observed site but no alternate allele calls (AC = 0)"
    } else if (af >= 0.05) {
      strength <- "strong"; q <- 1.0; dirn <- -1L
      txt <- "common (AF >= 5%; meets the ACMG BA1 stand-alone benign frequency threshold for most Mendelian disorders, though not for risk alleles)"
    } else if (af >= 0.01) {
      strength <- "strong"; q <- 0.9; dirn <- -1L
      txt <- "low-frequency (1-5%); generally too frequent for a highly penetrant rare Mendelian allele"
    } else if (af >= 1e-4) {
      strength <- "strong"; q <- 0.8; dirn <- 0L
      txt <- "uncommon (0.01-1%); interpretation depends on disease prevalence and penetrance"
    } else {
      strength <- "indirect"; q <- 0.8; dirn <- 1L
      txt <- "rare (< 0.01%); consistent with, but not proof of, pathogenicity"
    }
    if (s$an < CONFIG$gnomad_min_an) {
      low_an <- TRUE; q <- q * 0.6
      if (strength == "strong") strength <- "indirect"
    }
    if (!pass) q <- q * 0.5

    claim <- sprintf("gnomAD %s: %s allele frequency %s (AC %d / AN %d; %d homozygotes%s). Variant is %s.%s",
                     CONFIG$gnomad_dataset, vid, formatC(af, format = "g", digits = 3),
                     as.integer(s$ac), as.integer(s$an), as.integer(s$homozygotes),
                     if (!is.na(s$max_pop)) sprintf("; highest population AF %s in '%s'",
                                                    formatC(s$max_pop_af, format = "g", digits = 3), s$max_pop) else "",
                     txt, if (!pass) sprintf(" Quality filters failed: %s.", paste(s$filters, collapse = ", ")) else "")
    items[[length(items) + 1]] <- new_item(
      src, claim, strength, q, replicated = s$an >= CONFIG$gnomad_min_an,
      study_ids = paste0("gnomad:", vid), direction = dirn,
      details = list(variant_id = vid, dataset = CONFIG$gnomad_dataset, ac = s$ac, an = s$an,
                     af = af, homozygotes = s$homozygotes, max_population = s$max_pop,
                     max_population_af = s$max_pop_af, filters = as.list(s$filters),
                     sequencing = as.list(s$sequencing)),
      refs = paste0("gnomAD:", vid))
    refs <- c(refs, list(new_ref(vid, "database_record", src, "gnomAD",
                                 url = sprintf("https://gnomad.broadinstitute.org/variant/%s?dataset=%s", vid, CONFIG$gnomad_dataset))))
  }

  if (!length(items)) {
    if (!is.null(any_error)) return(module_result(src, "error", error = any_error))
    return(module_result(src, "no_data"))
  }
  if (!is.null(any_error)) notes <- c(notes, paste("Partial gnomAD failure:", any_error))
  module_result(src, "ok", items, refs, notes,
                flags = list(low_an = low_an, filter_fail = filter_fail),
                raw = list(max_af = max_af))
}

# -----------------------------------------------------------------------------
# 3. GWAS Catalog (REST v1 with v2 fallback)
# -----------------------------------------------------------------------------
parse_sample_size <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_real_)
  nums <- regmatches(s, gregexpr("[0-9][0-9,]*", s))[[1]]
  nums <- suppressWarnings(as.numeric(gsub(",", "", nums)))
  if (!length(nums) || all(is.na(nums))) NA_real_ else sum(nums, na.rm = TRUE)
}

extract_traits <- function(x) {
  if (is.null(x) || !length(x)) return(character())
  out <- vapply(x, function(t) {
    if (is.list(t)) chr1(get_any(t, c("trait", "efo_trait", "reported_trait", "label")), "")
    else as.character(t)
  }, "")
  out[nzchar(out)]
}

effect_sign <- function(or, beta, beta_dir) {
  if (!is.na(or) && or > 0 && or != 1) return(sign(log(or)))
  if (!is.na(beta) && beta != 0) {
    d <- tolower(if (is.na(beta_dir)) "" else beta_dir)
    if (grepl("decrease", d)) return(-1)
    if (grepl("increase", d)) return(1)
    return(sign(beta))
  }
  NA_real_
}

gwas_study_lookup <- function(a, api, cache) {
  href <- chr1(pluck(a, "_links", "study", "href"))
  acc <- chr1(get_any(a, c("accession_id", "accessionId", "study_accession")))
  key <- if (!is.na(href)) href else acc
  if (is.na(key)) return(NULL)
  if (!is.null(cache$studies[[key]])) return(cache$studies[[key]])
  if (cache$n >= CONFIG$gwas_max_study_lookups) { cache$truncated <- TRUE; return(NULL) }
  cache$n <- cache$n + 1
  url <- if (!is.na(href)) href else sprintf("%s/studies/%s", CONFIG$gwas_v2_api, acc)
  r <- http_request(url, realm = "ebi", rate = 5)
  st <- if (r$ok) r$data else list()
  cache$studies[[key]] <- st
  st
}

normalize_gwas_assoc <- function(a, api, cache) {
  study <- a$study
  if (is.null(study)) study <- gwas_study_lookup(a, api, cache) %||% list()
  pub <- study$publicationInfo %||% list()

  p <- num1(get_any(a, c("pvalue", "p_value")))
  if (is.na(p)) {
    m <- num1(get_any(a, c("pvalueMantissa", "pvalue_mantissa")))
    e <- num1(get_any(a, c("pvalueExponent", "pvalue_exponent")))
    if (!is.na(m) && !is.na(e)) p <- m * 10^e
  }
  risk <- character()
  for (l in a$loci %||% list()) for (ra in l$strongestRiskAlleles %||% list())
    risk <- c(risk, chr1(ra$riskAlleleName, ""))
  if (!length(risk)) risk <- as.character(unlist(get_any(a, c("snp_effect_allele", "risk_allele",
                                                                "strongest_risk_allele"))))
  allele <- toupper(sub("^.*-", "", risk[1]))

  traits <- extract_traits(get_any(a, c("efoTraits", "efo_traits")))
  if (!length(traits)) traits <- chr1(get_any(a, c("reported_trait", "traitName")),
                                      chr1(pluck(study, "diseaseTrait", "trait"), "unspecified trait"))
  init_ss <- chr1(get_any(study, c("initialSampleSize", "initial_sample_size")), "")
  rep_ss <- chr1(get_any(study, c("replicationSampleSize", "replication_sample_size")), "")

  list(
    trait = paste(traits, collapse = "; "),
    pvalue = p,
    allele = allele,
    or = num1(get_any(a, c("orPerCopyNum", "or_value", "or_per_copy_number"))),
    beta = num1(get_any(a, c("betaNum", "beta"))),
    beta_dir = chr1(get_any(a, c("betaDirection", "beta_direction"))),
    study_acc = chr1(get_any(study, c("accessionId", "accession_id")),
                     chr1(get_any(a, c("accession_id", "accessionId")))),
    pmid = chr1(get_any(pub, c("pubmedId")),
                chr1(get_any(study, c("pubmed_id")), chr1(get_any(a, c("pubmed_id", "pubmedId"))))),
    title = chr1(get_any(pub, c("title")), chr1(get_any(study, c("title")))),
    year = year_of(get_any(pub, c("publicationDate")) %||% get_any(study, c("publication_date"))),
    journal = chr1(get_any(pub, c("publication")), chr1(get_any(study, c("journal")))),
    initial_ss = init_ss,
    n_initial = parse_sample_size(init_ss),
    has_replication = nzchar(trimws(rep_ss)) &&
      !grepl("^(NA|N/A|none|-)$", trimws(rep_ss), ignore.case = TRUE) &&
      isTRUE(parse_sample_size(rep_ss) > 0)
  )
}

query_gwas <- function(var) {
  src <- "GWAS Catalog"
  rsid <- var$rsid_current
  r1 <- http_request(sprintf("%s/singleNucleotidePolymorphisms/%s/associations", CONFIG$gwas_v1_api, rsid),
                     query = list(projection = "associationBySnp"), realm = "ebi", rate = 5)
  api <- "v1"; assoc <- if (r1$ok) pluck(r1$data, "_embedded", "associations") %||% list() else list()
  if (!length(assoc)) {
    r2 <- http_request(paste0(CONFIG$gwas_v2_api, "/associations"),
                       query = list(rs_id = rsid, size = 500), realm = "ebi", rate = 5)
    if (r2$ok) { api <- "v2"; assoc <- pluck(r2$data, "_embedded", "associations") %||% list() }
    if (!length(assoc)) {
      if (r1$ok || r2$ok || isTRUE(r1$status == 404)) {
        return(module_result(src, "no_data", notes = sprintf("No GWAS Catalog associations for %s.", rsid)))
      }
      return(module_result(src, "error", error = r1$error %||% r2$error))
    }
  }

  cache <- new.env(); cache$studies <- list(); cache$n <- 0; cache$truncated <- FALSE
  rows <- lapply(assoc, normalize_gwas_assoc, api = api, cache = cache)
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  df$trait_key <- tolower(trimws(df$trait))
  df$sign <- mapply(effect_sign, df$or, df$beta, df$beta_dir)
  # orient effects to the first alternate allele
  df$sign[!is.na(df$allele) & df$allele == var$ref] <- -df$sign[!is.na(df$allele) & df$allele == var$ref]
  df$sign[is.na(df$allele) | !(df$allele %in% c(var$ref, var$alts))] <- NA

  ss <- tolower(unique(df$initial_ss[nzchar(df$initial_ss)]))
  non_eur <- "african|asian|hispanic|latin|admixed|japanese|chinese|korean|indian|arab|hawaiian|native|oceanian|middle eastern|mexican|brazilian|multi-ancestry|trans-ancestry"
  eur_only <- length(ss) > 0 && any(grepl("european", ss)) && !any(grepl(non_eur, ss))

  groups <- split(df, df$trait_key)
  best_p <- vapply(groups, function(g) safe_max(-log10(g$pvalue)), numeric(1))
  groups <- groups[order(-best_p)]
  n_traits_total <- length(groups)
  groups <- head(groups, CONFIG$gwas_max_trait_items)

  gwas_summ <- function(rows) {
    sig <- rows[!is.na(rows$pvalue) & rows$pvalue < CONFIG$gwas_sig_p, , drop = FALSE]
    use <- if (nrow(sig)) sig else rows
    keys <- unique(ifelse(!is.na(use$pmid), paste0("pmid:", use$pmid), paste0("gcst:", use$study_acc)))
    keys <- keys[!grepl(":NA$", keys)]
    list(sig = nrow(sig) > 0, use = use, study_keys = keys, n_ind = length(keys),
         replicated = length(keys) >= 2 || any(use$has_replication),
         meta = any(grepl("meta-analys", tolower(use$title))),
         max_n = safe_max(use$n_initial),
         min_p = suppressWarnings(min(use$pvalue, na.rm = TRUE)),
         signs = unique(stats::na.omit(use$sign)))
  }
  gwas_assertions <- function(rows, trait) {
    lapply(seq_len(nrow(rows)), function(k) {
      sg <- rows$sign[k]
      is_meta <- grepl("meta-analys", tolower(chr1(rows$title[k], "")))
      who <- if (!is.na(rows$pmid[k])) paste0("PMID:", rows$pmid[k]) else chr1(rows$study_acc[k], paste0("row", k))
      new_assertion(
        id = sprintf("%s#%d", who, k), source = "GWAS Catalog", group = trait,
        category = if (is.na(sg)) NA_character_ else if (sg > 0) "positive_effect" else "negative_effect",
        date = if (is.na(rows$year[k])) as.Date(NA) else as.Date(sprintf("%d-07-01", as.integer(rows$year[k]))),
        tier = if (is_meta) 3 else if (isTRUE(rows$has_replication[k])) 2 else 1,
        submitter = who, n = rows$n_initial[k], meta = is_meta,
        label = sprintf("%s (%s; p = %s)", chr1(rows$title[k], who), chr1(rows$year[k], "n.d."),
                        fmt_p(rows$pvalue[k])))
    })
  }

  items <- list(); refs <- list(); small_n <- FALSE; conflicts <- list()
  for (gk in names(groups)) {
    m <- gwas_summ(groups[[gk]])
    all_rows <- m$use
    trait <- all_rows$trait[1]
    conflict <- length(m$signs) > 1
    ctx <- "baseline_no_conflict"; superseded_note <- ""

    if (conflict) {
      cf <- resolve_conflict(gwas_assertions(all_rows, trait), src, trait)
      if (identical(cf$status, "resolved_correction")) {
        keep_sign <- if (cf$current_position == "positive_effect") 1 else -1
        m <- gwas_summ(all_rows[is.na(all_rows$sign) | all_rows$sign == keep_sign, , drop = FALSE])
        conflict <- FALSE
        ctx <- "corrected_by_review"
        superseded_note <- sprintf("; earlier opposite-direction report(s) (%s) superseded by later evidence (%s review: %s)",
                                   paste(head(vapply(cf$superseded, `[[`, "", "id"), 4), collapse = ", "),
                                   cf$review_signal, paste(unlist(cf$signal_basis), collapse = "; "))
        conflicts[[length(conflicts) + 1]] <- cf
      } else if (identical(cf$status, "unresolved")) {
        ctx <- "unresolved_conflict"
        conflicts[[length(conflicts) + 1]] <- cf
      }
    }
    use <- m$use

    if (conflict) {
      strength <- "mixed"; q <- 1.0
    } else if (m$sig) {
      strength <- "direct"
      q <- 1.0 + 0.2 * m$replicated + 0.2 * m$meta
      if (!is.na(m$max_n) && m$max_n >= CONFIG$gwas_large_n) q <- q + 0.1
      if (!is.na(m$max_n) && m$max_n < CONFIG$gwas_small_n) { q <- q - 0.3; small_n <- TRUE }
    } else {
      strength <- "emerging"; q <- 0.6
    }

    claim <- sprintf("%s association with '%s' (best p = %s; %d association record(s) from %d independent publication(s)%s%s%s%s%s).",
                     if (m$sig) "Genome-wide significant" else "Sub-threshold (suggestive)",
                     trait, fmt_p(m$min_p), nrow(use), m$n_ind,
                     if (m$replicated) "; replicated" else "; not replicated",
                     if (m$meta) "; includes meta-analysis" else "",
                     if (!is.na(m$max_n)) sprintf("; largest discovery N ~ %s", format(m$max_n, big.mark = ",")) else "",
                     if (conflict) "; effect directions disagree across studies with no reviewed resolution" else "",
                     superseded_note)

    pm <- unique(stats::na.omit(all_rows$pmid))
    acc <- unique(stats::na.omit(all_rows$study_acc))
    for (p in pm) {
      row <- all_rows[!is.na(all_rows$pmid) & all_rows$pmid == p, , drop = FALSE][1, ]
      refs <- c(refs, list(new_ref(p, "pmid", src, claimed = list(title = row$title, year = row$year,
                                                                   journal = row$journal))))
    }
    for (a in acc) refs <- c(refs, list(new_ref(a, "database_record", src, "GWAS Catalog",
                                                url = sprintf("https://www.ebi.ac.uk/gwas/studies/%s", a))))
    items[[length(items) + 1]] <- new_item(
      src, claim, strength, q, m$replicated, study_ids = m$study_keys,
      direction = if (m$sig && !conflict) 1L else 0L,
      details = list(trait = trait, best_p = m$min_p, n_records = nrow(use),
                     n_independent_publications = m$n_ind, replicated = m$replicated,
                     meta_analysis = m$meta, max_discovery_n = m$max_n,
                     effect_direction_conflict = ctx != "baseline_no_conflict",
                     conflict_status = ctx, api_version = api),
      refs = c(paste0("PMID:", pm), paste0("GWAS Catalog:", acc)),
      conflict_context = ctx,
      conflict_key = if (ctx == "baseline_no_conflict") NA_character_ else paste(src, trait, sep = "|"))
  }

  notes <- if (n_traits_total > length(groups))
    sprintf("%d additional traits not itemised (only the top %d by p-value are scored).",
            n_traits_total - length(groups), length(groups)) else character()
  module_result(src, "ok", items, refs, notes,
                flags = list(eur_only = eur_only, small_n = small_n,
                             truncated = n_traits_total > length(groups) || isTRUE(cache$truncated)),
                raw = list(n_associations = nrow(df), n_traits = n_traits_total, api = api,
                           conflicts = conflicts))
}

# -----------------------------------------------------------------------------
# 4. AlphaFold DB (+ AlphaMissense annotations)
# -----------------------------------------------------------------------------
residue_plddt <- function(pdb_text, resnum) {
  lines <- strsplit(pdb_text, "\n", fixed = TRUE)[[1]]
  ca <- lines[startsWith(lines, "ATOM") & substr(lines, 13, 16) == " CA "]
  if (!length(ca)) return(NA_real_)
  rn <- suppressWarnings(as.integer(substr(ca, 23, 26)))
  bf <- suppressWarnings(as.numeric(substr(ca, 61, 66)))
  v <- bf[!is.na(rn) & rn == resnum]
  if (length(v)) v[1] else NA_real_
}

am_direction <- function(cls, score) {
  c <- if (is.na(cls)) "" else tolower(cls)
  if (grepl("path", c)) return(1L)
  if (grepl("ben", c)) return(-1L)
  if (nzchar(c)) return(0L)
  if (!is.na(score)) return(if (score > 0.564) 1L else if (score < 0.34) -1L else 0L)
  0L
}

query_alphafold <- function(var) {
  src <- "AlphaFold"
  if (is.na(var$uniprot)) {
    return(module_result(src, "skipped", notes = "No Swiss-Prot accession on the selected transcript (likely non-coding); AlphaFold/AlphaMissense assess protein positions only."))
  }
  if (is.na(var$protein_pos)) {
    return(module_result(src, "skipped", notes = "Variant does not map to a protein position on the selected transcript."))
  }
  r <- http_request(sprintf("%s/prediction/%s", CONFIG$alphafold_api, var$uniprot), realm = "afdb", rate = 5)
  if (!r$ok) {
    if (isTRUE(r$status == 404)) return(module_result(src, "no_data", notes = sprintf("No AlphaFold model for %s.", var$uniprot)))
    return(module_result(src, "error", error = r$error))
  }
  entries <- r$data
  if (!is.null(names(entries))) entries <- list(entries)
  pos <- var$protein_pos
  ent <- NULL
  for (e in entries) {
    if (!identical(chr1(e$uniprotAccession), var$uniprot)) next
    if (pos >= num1(e$uniprotStart, 1) && pos <= num1(e$uniprotEnd, Inf)) { ent <- e; break }
  }
  if (is.null(ent)) ent <- entries[[1]]

  entry_id <- chr1(ent$entryId, paste0("AF-", var$uniprot, "-F1"))
  model_version <- chr1(ent$latestVersion, chr1(ent$modelCreatedDate, "unknown"))
  global_plddt <- num1(ent$globalMetricValue)
  local_pos <- as.integer(pos - num1(ent$uniprotStart, 1) + 1)

  useq <- chr1(ent$uniprotSequence, "")
  seq_aa <- if (nzchar(useq) && local_pos >= 1 && local_pos <= nchar(useq)) substr(useq, local_pos, local_pos) else NA_character_
  mismatch <- !is.na(seq_aa) && !is.na(var$ref_aa) && seq_aa != var$ref_aa
  db_flags <- if (mismatch) "reference_residue_mismatch" else character()
  refs <- list(new_ref(entry_id, "database_record", src, "AlphaFold DB", flags = db_flags,
                       url = sprintf("https://alphafold.ebi.ac.uk/entry/%s", var$uniprot)))
  ref_tag <- paste0("AlphaFold DB:", entry_id)
  items <- list(); notes <- character()
  if (mismatch) notes <- c(notes, sprintf("Transcript reference residue (%s) differs from UniProt residue %d (%s): possible isoform mismatch; AlphaFold items down-weighted.", var$ref_aa, as.integer(pos), seq_aa))

  plddt <- NA_real_
  pdb_url <- chr1(ent$pdbUrl)
  if (!is.na(pdb_url)) {
    p <- http_request(pdb_url, parse = "text", realm = "afdb", rate = 5)
    if (p$ok) plddt <- residue_plddt(p$data, local_pos)
  }
  if (!is.na(plddt)) {
    if (plddt >= 90)      { strength <- "indirect"; q <- 1.0; lab <- "very high" }
    else if (plddt >= 70) { strength <- "indirect"; q <- 0.8; lab <- "confident" }
    else if (plddt >= 50) { strength <- "emerging"; q <- 0.6; lab <- "low" }
    else                  { strength <- "emerging"; q <- 0.4; lab <- "very low (likely disordered)" }
    if (mismatch) q <- q * 0.5
    items[[length(items) + 1]] <- new_item(
      src, sprintf("Residue %s%d of %s (%s) lies in a %s-confidence region of AlphaFold model %s (residue pLDDT %.1f; model mean %.1f). Predicted structure; not experimentally determined.",
                   chr1(var$ref_aa, "?"), as.integer(pos), var$uniprot, chr1(var$gene, "?"), lab,
                   entry_id, plddt, global_plddt),
      strength, q, study_ids = paste0("afdb:", entry_id), direction = 0L,
      details = list(entry_id = entry_id, model_version = model_version, residue = as.integer(pos),
                     residue_plddt = plddt, global_plddt = global_plddt,
                     reference_residue_check = if (is.na(seq_aa)) "not_checked" else if (mismatch) "mismatch" else "match"),
      refs = ref_tag)
  }

  am_res <- list()
  am_url <- chr1(ent$amAnnotationsUrl)
  if (length(var$missense) && !is.na(am_url)) {
    a <- http_request(am_url, parse = "text", realm = "afdb", rate = 5)
    am <- if (a$ok) tryCatch(utils::read.csv(text = a$data, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
    if (!is.null(am) && all(c("protein_variant", "am_pathogenicity") %in% names(am))) {
      for (ms in var$missense) {
        hit <- am[am$protein_variant == ms$protein_change, , drop = FALSE]
        if (!nrow(hit)) next
        score <- suppressWarnings(as.numeric(hit$am_pathogenicity[1]))
        cls <- if ("am_class" %in% names(hit)) as.character(hit$am_class[1]) else NA_character_
        dirn <- am_direction(cls, score)
        q <- if (!is.na(plddt) && plddt < 70) 0.6 else 1.0
        if (mismatch) q <- q * 0.5
        am_res[[length(am_res) + 1]] <- list(change = ms$protein_change, score = score, class = cls, direction = dirn)
        items[[length(items) + 1]] <- new_item(
          src, sprintf("AlphaMissense predicts p.%s as %s (pathogenicity score %.3f). Unsupervised computational prediction; hypothesis-generating only.",
                       ms$protein_change, if (is.na(cls)) c("likely benign", "ambiguous", "likely pathogenic")[dirn + 2] else cls, score),
          "emerging", q, study_ids = paste0("alphamissense:", ms$protein_change), direction = dirn,
          details = list(protein_change = ms$protein_change, am_pathogenicity = score, am_class = cls),
          refs = ref_tag)
      }
    } else {
      notes <- c(notes, "AlphaMissense annotations could not be retrieved or parsed.")
    }
  }
  if (!length(items)) return(module_result(src, "no_data", refs = list(), notes = c(notes, "No residue-level AlphaFold data retrieved.")))
  module_result(src, "ok", items, refs, notes, flags = list(isoform_mismatch = mismatch),
                raw = list(model_version = model_version, alphamissense = am_res, plddt = plddt))
}

# -----------------------------------------------------------------------------
# 5. JASPAR (PFM download + ref/alt PWM scan around the SNP)
# -----------------------------------------------------------------------------
pfm_to_pwm <- function(m, pseudo = 0.8) {
  p <- sweep(m + pseudo / 4, 2, colSums(m) + pseudo, "/")
  log2(p / 0.25)
}

parse_jaspar_pfms <- function(lines) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  hdr <- which(startsWith(lines, ">"))
  motifs <- vector("list", length(hdr))
  for (i in seq_along(hdr)) {
    h <- strsplit(sub("^>", "", lines[hdr[i]]), "\\s+")[[1]]
    block <- lines[hdr[i] + 1:4]
    if (anyNA(block)) next
    rows <- lapply(block, function(r) {
      nums <- gsub("\\[|\\]", " ", sub("^[ACGT]", "", r))
      suppressWarnings(as.numeric(strsplit(trimws(nums), "\\s+")[[1]]))
    })
    if (length(unique(lengths(rows))) != 1 || anyNA(unlist(rows))) next
    m <- do.call(rbind, rows)
    rownames(m) <- substr(block, 1, 1)
    if (!setequal(rownames(m), c("A", "C", "G", "T"))) next
    m <- m[c("A", "C", "G", "T"), , drop = FALSE]
    pwm <- pfm_to_pwm(m)
    motifs[[i]] <- list(id = h[1], name = if (length(h) > 1) paste(h[-1], collapse = " ") else h[1],
                        pwm = pwm, width = ncol(pwm),
                        min = sum(apply(pwm, 2, min)), max = sum(apply(pwm, 2, max)))
  }
  Filter(Negate(is.null), motifs)
}

load_jaspar_motifs <- function() {
  dir.create(CONFIG$cache_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(CONFIG$cache_dir, basename(CONFIG$jaspar_pfm_url))
  if (!file.exists(f)) {
    r <- http_request(CONFIG$jaspar_pfm_url, parse = "text", realm = "jaspar", rate = 2)
    if (!r$ok || !grepl(">MA", r$data, fixed = TRUE)) return(NULL)
    writeLines(r$data, f)
  }
  parse_jaspar_pfms(readLines(f, warn = FALSE))
}

encode_dna <- function(s) match(strsplit(toupper(s), "")[[1]], c("A", "C", "G", "T"))

motif_best_rel <- function(codes, center, mot) {
  w <- mot$width; L <- length(codes)
  lo <- max(1L, center - w + 1L); hi <- min(center, L - w + 1L)
  if (lo > hi) return(NA_real_)
  cols <- seq_len(w)
  sc <- vapply(lo:hi, function(s) {
    idx <- codes[s:(s + w - 1L)]
    if (anyNA(idx)) NA_real_ else sum(mot$pwm[cbind(idx, cols)])
  }, numeric(1))
  if (all(is.na(sc))) return(NA_real_)
  (max(sc, na.rm = TRUE) - mot$min) / (mot$max - mot$min)
}

scan_allele <- function(seq_chr, center, motifs) {
  f <- encode_dna(seq_chr)
  rc <- rev(5L - f)                      # A(1)<->T(4), C(2)<->G(3); N stays NA
  rc_center <- length(f) - center + 1L
  vapply(motifs, function(m) {
    s <- c(motif_best_rel(f, center, m), motif_best_rel(rc, rc_center, m))
    safe_max(s)
  }, numeric(1))
}

query_jaspar <- function(var) {
  src <- "JASPAR"
  if (!isTRUE(var$is_snv)) {
    return(module_result(src, "skipped", notes = "Motif scanning is implemented for single-nucleotide variants only."))
  }
  motifs <- load_jaspar_motifs()
  if (is.null(motifs) || !length(motifs)) {
    return(module_result(src, "error", error = "Could not download or parse JASPAR PFMs."))
  }
  flank <- as.integer(max(CONFIG$jaspar_flank, max(vapply(motifs, `[[`, numeric(1), "width"))))
  r <- ensembl_get(sprintf("/sequence/region/human/%s:%d..%d:1", var$chr,
                           as.integer(var$start - flank), as.integer(var$end + flank)))
  seq <- toupper(chr1(pluck(r$data, "seq"), ""))
  if (!nzchar(seq)) return(module_result(src, "error", error = "Could not retrieve flanking sequence from Ensembl."))
  center <- flank + 1L
  ref_obs <- substr(seq, center, center)
  ref_ok <- identical(ref_obs, var$ref)
  notes <- character()
  if (!ref_ok) notes <- c(notes, sprintf("Genome reference base (%s) does not match the dbSNP reference allele (%s); motif results down-weighted.", ref_obs, var$ref))

  ref_scores <- scan_allele(seq, center, motifs)
  hits <- list()
  for (alt in var$alts) {
    alt_seq <- paste0(substr(seq, 1, center - 1), alt, substr(seq, center + 1, nchar(seq)))
    alt_scores <- scan_allele(alt_seq, center, motifs)
    delta <- alt_scores - ref_scores
    keep <- which(!is.na(delta) & pmax(ref_scores, alt_scores) >= CONFIG$jaspar_min_rel &
                    abs(delta) >= CONFIG$jaspar_min_delta)
    for (k in keep) hits[[length(hits) + 1]] <- list(
      matrix_id = motifs[[k]]$id, tf = motifs[[k]]$name, alt = alt,
      ref_rel = ref_scores[k], alt_rel = alt_scores[k], delta = delta[k])
  }

  if (!length(hits)) {
    return(module_result(src, "no_data", notes = c(notes, sprintf(
      "No %s CORE vertebrate motif (of %d scanned) gained or lost a match (relative score >= %.2f, |delta| >= %.2f).",
      CONFIG$jaspar_release, length(motifs), CONFIG$jaspar_min_rel, CONFIG$jaspar_min_delta))))
  }
  hits <- hits[order(-abs(vapply(hits, `[[`, numeric(1), "delta")))]
  n_hits <- length(hits)
  hits <- head(hits, CONFIG$jaspar_max_items)

  items <- list(); refs <- list()
  db_flags <- if (!ref_ok) "reference_base_mismatch" else character()
  for (h in hits) {
    q <- 0.5 + 0.3 * (max(h$ref_rel, h$alt_rel) >= 0.9) + 0.2 * (var$regulatory_overlaps > 0)
    if (!ref_ok) q <- q * 0.5
    items[[length(items) + 1]] <- new_item(
      src, sprintf("In-silico %s of a %s motif (%s, %s): relative PWM score %.2f (ref %s) -> %.2f (alt %s). Sequence-only prediction without chromatin or tissue context%s.",
                   if (h$delta < 0) "loss" else "gain", h$tf, h$matrix_id, CONFIG$jaspar_release,
                   h$ref_rel, var$ref, h$alt_rel, h$alt,
                   if (var$regulatory_overlaps > 0) "; variant overlaps an Ensembl regulatory feature" else ""),
      "emerging", q, study_ids = paste0("jaspar:", h$matrix_id), direction = 0L,
      details = list(matrix_id = h$matrix_id, tf = h$tf, alt_allele = h$alt,
                     ref_relative_score = round(h$ref_rel, 3), alt_relative_score = round(h$alt_rel, 3),
                     delta = round(h$delta, 3), effect = if (h$delta < 0) "predicted_loss" else "predicted_gain"),
      refs = paste0("JASPAR:", h$matrix_id))
    refs <- c(refs, list(new_ref(h$matrix_id, "database_record", src, "JASPAR", flags = db_flags,
                                 url = sprintf("https://jaspar.elixir.no/matrix/%s/", h$matrix_id))))
  }
  if (n_hits > length(hits)) notes <- c(notes, sprintf("%d motif changes found; top %d by |delta| itemised.", n_hits, length(hits)))
  module_result(src, "ok", items, refs, notes, flags = list(ref_mismatch = !ref_ok),
                raw = list(n_motifs_scanned = length(motifs), n_hits = n_hits))
}

# -----------------------------------------------------------------------------
# 6. NIH Genetic Testing Registry (gene-level)
# -----------------------------------------------------------------------------
query_gtr <- function(var) {
  src <- "NIH GTR"
  if (is.na(var$gene)) return(module_result(src, "skipped", notes = "No gene symbol available for GTR lookup."))
  s <- eutils("esearch", list(db = "gtr", term = paste0(var$gene, "[gene]"), retmax = CONFIG$gtr_max_records))
  if (!s$ok) return(module_result(src, "error", error = s$error))
  count <- num1(pluck(s$data, "esearchresult", "count"), 0)
  ids <- as.character(unlist(pluck(s$data, "esearchresult", "idlist")))
  if (!length(ids)) return(module_result(src, "no_data", notes = sprintf("No GTR tests registered for %s.", var$gene)))

  recs <- list()
  sm <- eutils("esummary", list(db = "gtr", id = paste(ids, collapse = ",")))
  if (sm$ok) {
    res <- pluck(sm$data, "result")
    uids <- as.character(unlist(res$uids))
    recs <- Filter(function(x) is.list(x) && is.null(x$error), lapply(uids, function(u) {
      x <- res[[u]]; if (is.list(x)) x$.uid <- u; x
    }))
  }
  test_type <- vapply(recs, function(r) tolower(chr1(get_any(r, c("testtype", "test_type")), "")), "")
  n_clin <- if (length(recs)) sum(grepl("clinical", test_type) | !nzchar(test_type)) else count
  conds <- unlist(lapply(recs, function(r) vapply(get_any(r, c("conditionlist", "conditions")) %||% list(),
                                                  function(c) if (is.list(c)) chr1(get_any(c, c("name", "conditionname")), "") else as.character(c), "")))
  conds <- conds[nzchar(conds)]
  top_conds <- if (length(conds)) names(sort(table(conds), decreasing = TRUE))[seq_len(min(5, length(unique(conds))))] else character()
  accs <- vapply(recs, function(r) chr1(get_any(r, c("accession", "testaccession", "accessionversion")), paste0("uid", r$.uid)), "")
  accs <- head(accs, 10)

  refs <- lapply(seq_along(accs), function(i) new_ref(accs[i], "database_record", src, "NIH GTR",
                                                       url = sprintf("https://www.ncbi.nlm.nih.gov/gtr/tests/%s/", recs[[i]]$.uid)))
  q <- if (n_clin >= 5) 0.8 else 0.5
  item <- new_item(
    src, sprintf("%d genetic test record(s) for %s are registered in GTR (%d clinical in the %d summarised)%s. Gene-level evidence of clinical testing practice; not specific to this variant.",
                 as.integer(count), var$gene, as.integer(n_clin), length(recs),
                 if (length(top_conds)) paste0("; most frequent conditions: ", paste(top_conds, collapse = "; ")) else ""),
    "indirect", q, study_ids = paste0("gtr:", var$gene), direction = 0L,
    details = list(gene = var$gene, n_tests = count, n_clinical_in_sample = n_clin,
                   top_conditions = as.list(top_conds)),
    refs = paste0("NIH GTR:", accs))
  module_result(src, "ok", list(item), refs, raw = list(n_tests = count))
}

# -----------------------------------------------------------------------------
# Cross-source consistency checks
# -----------------------------------------------------------------------------
cross_source_checks <- function(mods) {
  items <- list(); conflicts <- list()
  recs <- mods$clinvar$raw$records %||% list()
  max_af <- mods$gnomad$raw$max_af %||% NA_real_
  am <- mods$alphafold$raw$alphamissense %||% list()

  # ClinVar records corrected to benign/VUS no longer have direction 1, so a
  # resolved correction removes this discordance automatically.
  plp <- Filter(function(r) r$direction == 1L && r$stars >= 1 && !r$conflict &&
                  !grepl("risk|association|protective|drug|affects", tolower(r$description)), recs)
  if (length(plp) && !is.na(max_af) && max_af >= 0.05) {
    key <- "ClinVar P/LP vs gnomAD population frequency"
    accs <- vapply(plp, `[[`, "", "accession")
    pd <- do.call(c, lapply(plp, function(r) safe_date(r$latest_plp_date)))
    predates <- length(pd) > 0 && !anyNA(pd) && all(pd < as.Date(CONFIG$gnomad_release_date))
    items[[length(items) + 1]] <- new_item(
      "Cross-source", sprintf("Discordance: ClinVar P/LP assertion(s) (%s) versus a gnomAD allele frequency of %s, too common for a highly penetrant Mendelian allele.%s Consider reduced penetrance, risk-allele status, or misclassification.",
                              paste(accs, collapse = ", "), formatC(max_af, format = "g", digits = 3),
                              if (predates) sprintf(" Every P/LP assertion was last evaluated before the gnomAD %s release (%s) and may predate current population data, but no newer reviewed reclassification exists, so the conflict is unresolved.",
                                                    CONFIG$gnomad_dataset, CONFIG$gnomad_release_date)
                              else " No newer reviewed reclassification resolves it."),
      "mixed", 1.0, direction = 0L, study_ids = "cross:clinvar_gnomad",
      refs = paste0("ClinVar:", accs), conflict_context = "unresolved_conflict",
      conflict_key = paste("Cross-source", key, sep = "|"))
    conflicts[[length(conflicts) + 1]] <- list(
      source = "Cross-source", group = key, status = "unresolved",
      positions = list("pathogenic/likely pathogenic (ClinVar)",
                       sprintf("common, max AF %s (gnomAD)", formatC(max_af, format = "g", digits = 3))),
      assertions_predate_population_data = predates,
      reason = "Population data contradict the clinical assertion, but no newer, reviewed ClinVar submission has reclassified it.")
  }

  ben <- Filter(function(r) r$direction == -1L && r$stars >= 2, recs)
  plp2 <- Filter(function(r) r$direction == 1L && r$stars >= 2, recs)
  am_p <- Filter(function(a) a$direction == 1L, am)
  am_b <- Filter(function(a) a$direction == -1L, am)
  if ((length(ben) && length(am_p)) || (length(plp2) && length(am_b))) {
    key <- "ClinVar (>= 2 stars) vs AlphaMissense"
    items[[length(items) + 1]] <- new_item(
      "Cross-source", "Discordance between a reviewed ClinVar classification (>= 2 stars) and the AlphaMissense prediction. Curated clinical evidence normally outweighs the computational prediction.",
      "mixed", 0.5, direction = 0L, study_ids = "cross:clinvar_alphamissense",
      conflict_context = "unresolved_conflict", conflict_key = paste("Cross-source", key, sep = "|"))
    conflicts[[length(conflicts) + 1]] <- list(
      source = "Cross-source", group = key, status = "unresolved",
      positions = list("ClinVar curated classification", "AlphaMissense prediction"),
      reason = "A computational prediction cannot review or correct a curated classification; recorded as an unresolved, low-weight disagreement.")
  }
  list(items = items, conflicts = conflicts)
}

# -----------------------------------------------------------------------------
# Confidence scoring
# -----------------------------------------------------------------------------
confidence_band <- function(score) {
  if (score >= 80) "High" else if (score >= 50) "Moderate" else if (score >= 25) "Low" else "Emerging"
}

# Booster points for conflicts resolved by a newer, reviewed position.
correction_boosters <- function(conflicts) {
  ids <- vapply(conflicts, function(cf) chr1(cf$id, ""), "")
  per <- stats::setNames(numeric(length(conflicts)), ids)
  res <- which(vapply(conflicts, function(cf) identical(cf$status, "resolved_correction"), logical(1)))
  if (!length(res)) return(list(total = 0, per = per))
  raw <- vapply(conflicts[res], function(cf)
    if (identical(cf$review_signal, "explicit")) CONFIG$booster_explicit else CONFIG$booster_implicit,
    numeric(1))
  o <- order(-raw)
  disc <- raw[o] * CONFIG$booster_decay^(seq_along(o) - 1)
  scale <- if (sum(disc) > CONFIG$booster_cap) CONFIG$booster_cap / sum(disc) else 1
  per[res[o]] <- disc * scale
  list(total = sum(per), per = per)
}

compute_confidence <- function(items, conflicts = list()) {
  pts <- CONFIG$strength_points
  method <- sprintf(paste(
    "Baseline: each item earns base points by strength (%s) x a quality multiplier (0-1.5).",
    "Within a source the i-th best item is discounted by %.1f^(i-1); total points T map to",
    "100*(1-exp(-T/%d)), multiplied by a consistency factor 1 - 0.5*(unresolved-conflict points / T).",
    "Concordant evidence receives baseline weight only. Booster: each conflict resolved by a newer,",
    "equal-or-higher-tier, reviewed position adds %d (explicit review) or %d (implicit review) points,",
    "discounted %.1f^(i-1) and capped at %d in total; superseded assertions carry no weight.",
    "Band gates are applied after the booster, so a correction cannot by itself lift a score past a gate."),
    paste(sprintf("%s=%d", names(pts), pts), collapse = ", "),
    CONFIG$within_source_decay, as.integer(CONFIG$saturation_k),
    as.integer(CONFIG$booster_explicit), as.integer(CONFIG$booster_implicit),
    CONFIG$booster_decay, as.integer(CONFIG$booster_cap))
  boost <- correction_boosters(conflicts)
  statuses <- vapply(conflicts, function(cf) chr1(cf$status, ""), "")
  n_res <- sum(statuses == "resolved_correction")
  n_unres <- sum(statuses == "unresolved")

  if (!length(items)) {
    return(list(score = 0L, band = "Emerging", band_definition = BAND_DEFS$Emerging,
                baseline_score = 0, correction_booster = 0,
                n_resolved_corrections = n_res, n_unresolved_conflicts = n_unres,
                raw_score = 0, consistency_factor = 1, total_weighted_points = 0,
                strength_counts = as.list(stats::setNames(rep(0L, length(pts)), names(pts))),
                gates_applied = list("No evidence items retrieved; score fixed at 0."),
                item_points = numeric(), booster_per_conflict = boost$per, method = method,
                rationale = "No evidence was retrieved from any source."))
  }
  src <- vapply(items, `[[`, "", "source")
  str <- vapply(items, `[[`, "", "strength")
  q <- vapply(items, `[[`, numeric(1), "quality")
  base <- unname(pts[str]) * q
  rank <- integer(length(items))
  for (s in unique(src)) {
    idx <- which(src == s)
    idx <- idx[order(-base[idx])]
    rank[idx] <- seq_along(idx)
  }
  eff <- base * CONFIG$within_source_decay^(rank - 1)
  total <- sum(eff)
  raw <- 100 * (1 - exp(-total / CONFIG$saturation_k))
  mixed_pts <- sum(eff[str == "mixed"])
  direct_pts <- sum(eff[str == "direct"])
  consistency <- if (total > 0) 1 - 0.5 * mixed_pts / total else 1
  baseline <- raw * consistency
  score <- min(100, baseline + boost$total)

  n <- table(factor(str, levels = names(pts)))
  rep_direct <- which(str == "direct" & vapply(items, function(i) isTRUE(i$replicated), logical(1)))
  indep <- unique(unlist(lapply(items[rep_direct], `[[`, "study_ids")))
  gates <- character()

  if (score >= 80 && !(length(rep_direct) >= 2 && length(indep) >= 2)) {
    score <- 79
    gates <- c(gates, "High requires >= 2 replicated direct items drawn from >= 2 independent studies/records; capped at 79.")
  }
  if (score >= 50 && !(n[["direct"]] >= 1 || n[["strong"]] >= 2)) {
    score <- 49
    gates <- c(gates, "Moderate requires >= 1 direct item or >= 2 strong items; capped at 49.")
  }
  only_weak <- (n[["direct"]] + n[["strong"]] + n[["indirect"]]) == 0
  conflict_dom <- mixed_pts > 0 && mixed_pts >= direct_pts
  if (score >= 25 && (only_weak || conflict_dom)) {
    score <- 24
    gates <- c(gates, if (only_weak) "Only emerging/mixed evidence available; capped at 24 (Emerging)."
               else "Unresolved conflicting evidence outweighs non-conflicting direct evidence; capped at 24 (Emerging).")
  }
  score <- as.integer(round(score))
  band <- confidence_band(score)

  rationale <- sprintf(
    "%d evidence item(s): %d direct (%d replicated), %d strong, %d indirect, %d mixed, %d emerging. Weighted points %.1f -> raw %.1f x consistency %.2f = baseline %.1f; correction booster +%.1f from %d resolved conflict(s) (%d unresolved) -> final %d (%s).%s",
    length(items), n[["direct"]], length(rep_direct), n[["strong"]], n[["indirect"]], n[["mixed"]],
    n[["emerging"]], total, raw, consistency, baseline, boost$total, n_res, n_unres, score, band,
    if (length(gates)) paste0(" Gates: ", paste(gates, collapse = " ")) else "")

  list(score = score, band = band, band_definition = BAND_DEFS[[band]],
       baseline_score = round(baseline, 1), correction_booster = round(boost$total, 1),
       n_resolved_corrections = n_res, n_unresolved_conflicts = n_unres,
       raw_score = round(raw, 1), consistency_factor = round(consistency, 3),
       total_weighted_points = round(total, 2),
       strength_counts = as.list(stats::setNames(as.integer(n), names(n))),
       independent_replicated_direct_studies = length(indep),
       gates_applied = as.list(gates), item_points = eff, booster_per_conflict = boost$per,
       method = method, rationale = rationale)
}

# -----------------------------------------------------------------------------
# Citation-integrity audit
# -----------------------------------------------------------------------------
title_similarity <- function(a, b) {
  stop_w <- c("", "the", "of", "and", "in", "a", "an", "for", "with", "to", "on", "by", "from", "at", "is")
  tok <- function(x) unique(setdiff(strsplit(gsub("[^a-z0-9 ]", " ", tolower(x)), "\\s+")[[1]], stop_w))
  A <- tok(a); B <- tok(b)
  if (!length(A) || !length(B)) return(NA_real_)
  length(intersect(A, B)) / length(union(A, B))
}

verify_pmids <- function(pmids) {
  out <- list()
  if (!length(pmids)) return(out)
  for (chunk in split(pmids, ceiling(seq_along(pmids) / 150))) {
    r <- eutils("esummary", list(db = "pubmed", id = paste(chunk, collapse = ",")))
    if (!r$ok) {
      for (p in chunk) out[[p]] <- list(lookup_ok = FALSE)
      next
    }
    res <- pluck(r$data, "result")
    for (p in chunk) {
      rec <- res[[p]]
      if (is.null(rec) || !is.null(rec$error)) { out[[p]] <- list(lookup_ok = TRUE, found = FALSE); next }
      ids <- rec$articleids %||% list()
      doi <- NA_character_
      for (x in ids) if (identical(chr1(x$idtype), "doi")) doi <- chr1(x$value)
      out[[p]] <- list(lookup_ok = TRUE, found = TRUE,
                       title = chr1(rec$title), journal = chr1(rec$fulljournalname, chr1(rec$source)),
                       year = year_of(rec$pubdate %||% rec$sortpubdate), doi = doi,
                       n_authors = length(rec$authors %||% list()),
                       pubtypes = as.character(unlist(rec$pubtype)))
    }
  }
  out
}

SEVERE_FLAGS <- c("invalid_identifier_format", "pmid_not_found_in_pubmed", "retracted_publication",
                  "is_retraction_notice", "expression_of_concern", "year_mismatch", "title_mismatch",
                  "reference_residue_mismatch", "reference_base_mismatch", "unexpected_accession_format")
# Informational flags are reported but do not downgrade integrity status
INFO_FLAGS <- c("no_doi_recorded", "no_source_metadata_to_cross_check")

audit_entry <- function(r, pm_info, ev_map) {
  flags <- r$flags
  retrieved <- list()
  resolved <- TRUE
  if (r$ref_type == "pmid") {
    if (!grepl("^[0-9]{1,9}$", r$raw_id)) {
      flags <- c(flags, "invalid_identifier_format"); resolved <- FALSE
    } else {
      info <- pm_info[[r$raw_id]]
      if (is.null(info) || !isTRUE(info$lookup_ok)) {
        flags <- c(flags, "pubmed_lookup_failed"); resolved <- NA
      } else if (!isTRUE(info$found)) {
        flags <- c(flags, "pmid_not_found_in_pubmed"); resolved <- FALSE
      } else {
        retrieved <- list(title = info$title, journal = info$journal, year = info$year,
                          doi = info$doi, n_authors = info$n_authors, publication_types = as.list(info$pubtypes))
        pt <- tolower(info$pubtypes)
        if ("retracted publication" %in% pt || grepl("^retracted", chr1(info$title, ""), ignore.case = TRUE))
          flags <- c(flags, "retracted_publication")
        if ("retraction of publication" %in% pt) flags <- c(flags, "is_retraction_notice")
        if (any(grepl("expression of concern", pt))) flags <- c(flags, "expression_of_concern")
        if (is.na(info$title) || is.na(info$journal) || is.na(info$year) || info$n_authors == 0)
          flags <- c(flags, "metadata_incomplete")
        if (is.na(info$doi)) flags <- c(flags, "no_doi_recorded")
        cl <- r$claimed
        cy <- suppressWarnings(as.integer(chr1(cl$year)))
        if (!is.na(cy) && !is.na(info$year) && abs(cy - info$year) > 1) flags <- c(flags, "year_mismatch")
        ct <- chr1(cl$title)
        if (!is.na(ct) && !is.na(info$title)) {
          sim <- title_similarity(ct, info$title)
          if (!is.na(sim) && sim < 0.5) flags <- c(flags, "title_mismatch")
        }
        if (is.na(ct) && is.na(cy)) flags <- c(flags, "no_source_metadata_to_cross_check")
      }
    }
  } else {
    pat <- ACC_PATTERNS[[r$database]]
    if (!is.null(pat) && !grepl(pat, r$raw_id)) flags <- c(flags, "unexpected_accession_format")
  }
  flags <- unique(flags)
  status <- if (is.na(resolved)) "unverified"
            else if (!resolved || any(flags %in% SEVERE_FLAGS)) "flagged"
            else if (length(setdiff(flags, INFO_FLAGS))) "verified_with_warnings"
            else "verified"
  list(ref_id = r$ref_id, ref_type = r$ref_type, database = r$database, url = r$url,
       cited_by = as.list(unique(r$cited_by)),
       evidence_ids = as.list(unique(ev_map[[r$ref_id]] %||% character())),
       integrity_status = status, flags = as.list(flags),
       check_performed = if (r$ref_type == "pmid")
         "PubMed esummary: existence, retraction/concern status, metadata completeness, title/year cross-check against citing source"
       else "Record retrieved live from the source API in this session; accession format validated",
       claimed_metadata = r$claimed, retrieved_metadata = retrieved)
}

build_audit <- function(items, refs) {
  if (!length(refs)) return(list())
  by_id <- list()
  for (r in refs) {
    cur <- by_id[[r$ref_id]]
    if (is.null(cur)) {
      by_id[[r$ref_id]] <- r
    } else {
      cur$cited_by <- union(cur$cited_by, r$cited_by)
      cur$flags <- union(cur$flags, r$flags)
      if (!length(cur$claimed) && length(r$claimed)) cur$claimed <- r$claimed
      by_id[[r$ref_id]] <- cur
    }
  }
  ev_map <- list()
  for (it in items) for (rid in it$refs) ev_map[[rid]] <- c(ev_map[[rid]], it$id)
  pm_refs <- Filter(function(r) r$ref_type == "pmid", by_id)
  pmids <- unique(vapply(pm_refs, `[[`, "", "raw_id"))
  pm_info <- verify_pmids(pmids[grepl("^[0-9]{1,9}$", pmids)])
  unname(lapply(by_id, audit_entry, pm_info = pm_info, ev_map = ev_map))
}

# -----------------------------------------------------------------------------
# Limitations (enumerated per query)
# -----------------------------------------------------------------------------
build_limitations <- function(var, mods, items, audit, retrieved_at, conflicts = list()) {
  L <- list()
  add <- function(code, description, sources = character()) {
    L[[length(L) + 1]] <<- list(id = sprintf("L%d", length(L) + 1), code = code,
                                description = description, affected_sources = as.list(sources))
  }
  st <- vapply(mods, `[[`, "", "status")
  srcs_of <- function(keys) unname(SOURCES[intersect(keys, names(mods))])
  item_src <- vapply(items, `[[`, "", "source")
  item_str <- vapply(items, `[[`, "", "strength")

  # knowledge cutoff / snapshot dates
  stale <- mods$clinvar$flags$stale %||% character()
  add("knowledge_cutoff", paste0(
    sprintf("Evidence reflects database contents retrieved at %s (gnomAD '%s', %s CORE vertebrate PFMs, AlphaFold DB model version %s; live ClinVar, GWAS Catalog and GTR records). ",
            retrieved_at, CONFIG$gnomad_dataset, CONFIG$jaspar_release,
            mods$alphafold$raw$model_version %||% "n/a"),
    "Curated resources lag the primary literature, so findings published after each resource's last curation cycle are not captured, and no full-text literature was read.",
    if (length(stale)) sprintf(" ClinVar classification(s) %s were last evaluated more than %d years ago.",
                               paste(stale, collapse = ", "), CONFIG$stale_years) else ""),
    srcs_of(names(mods)))

  # human clinical data
  if (!any(item_str == "direct")) {
    add("lack_of_human_clinical_data",
        "No direct human clinical or association evidence (a ClinVar assertion with review criteria, or a genome-wide-significant GWAS association) was found; remaining evidence is population-frequency, gene-level, or computational.",
        srcs_of(c("clinvar", "gwas")))
  }

  # sample sizes
  ss <- character()
  if (isTRUE(mods$gwas$flags$small_n)) ss <- c(ss, sprintf("at least one GWAS discovery cohort had N < %s", format(CONFIG$gwas_small_n, big.mark = ",")))
  if (isTRUE(mods$gnomad$flags$low_an)) ss <- c(ss, sprintf("gnomAD allele number < %s (low coverage or callability)", format(CONFIG$gnomad_min_an, big.mark = ",")))
  if (isTRUE(mods$clinvar$flags$single_submitter)) ss <- c(ss, "ClinVar assertion(s) rest on a single submitter")
  if (length(ss)) add("small_sample_sizes", paste0("Limited sample size or observations: ", paste(ss, collapse = "; "), "."),
                      srcs_of(c("gwas", "gnomad", "clinvar")))

  # RCTs
  add("absence_of_rcts",
      "None of the queried resources record randomised controlled trials. 'Direct' evidence here is observational (GWAS) or curated clinical classification, so the RCT/meta-analysis standard for 'High' is approximated by replicated or meta-analysed GWAS and multi-submitter/expert-panel ClinVar consensus.",
      srcs_of(names(mods)))

  # species / tissue
  if (any(item_src %in% c("JASPAR", "AlphaFold"))) {
    add("species_tissue_specificity", paste(
      "Computational items lack biological context:",
      if (any(item_src == "JASPAR")) "JASPAR profiles derive from assays (SELEX, ChIP-seq) across varied species and cell types, and predicted motif changes do not account for chromatin accessibility, TF expression, or tissue specificity." else "",
      if (any(item_src == "AlphaFold")) "AlphaFold models are single-chain predictions without ligands, partners, or post-translational modifications; AlphaMissense is a computational predictor without experimental validation." else ""),
      srcs_of(c("jaspar", "alphafold")))
  }

  # ancestry representation
  if (isTRUE(mods$gwas$flags$eur_only)) {
    add("population_ancestry_bias",
        "All GWAS discovery samples with stated ancestry are European; effect sizes and associations may not transfer to other ancestries.",
        srcs_of("gwas"))
  }

  # publication bias
  if (any(item_src %in% c("ClinVar", "GWAS Catalog"))) {
    add("publication_bias",
        "The GWAS Catalog only curates published associations reaching reporting thresholds, and ClinVar submissions skew toward reportable or clinically tested variants; null results are under-represented, which can inflate apparent consistency.",
        srcs_of(c("gwas", "clinvar")))
  }

  # conflict / temporal inference
  if (any(names(mods) %in% c("clinvar", "gwas"))) {
    st_c <- vapply(conflicts, function(cf) chr1(cf$status, ""), "")
    add("temporal_inference", paste0(
      sprintf("%d conflict(s) detected: %d resolved by a newer reviewed position, %d unresolved. ",
              length(st_c), sum(st_c == "resolved_correction"), sum(st_c == "unresolved")),
      sprintf("A correction is inferred only when a position postdates every opposing assertion by >= %d days, comes from an equal-or-higher review tier, and carries a review signal; newer alone never counts. ",
              as.integer(CONFIG$correction_min_gap_days)),
      "Dates are last-evaluated or publication dates, not the date a correction was made. ClinVar's API exposes only each submitter's current classification, so a laboratory that revised its own call is invisible and reads as agreement. ",
      "The GWAS Catalog omits null replications, so corrections by non-replication cannot be detected. ",
      "Concordant evidence is scored at baseline because agreement may reflect uncritical propagation; the tool cannot tell independent agreement from copying, so genuinely independent replication may be under-credited."),
      srcs_of(c("clinvar", "gwas")))
  }

  # citation metadata
  if (length(audit)) {
    stat <- vapply(audit, `[[`, "", "integrity_status")
    bad <- sum(stat != "verified")
    if (bad > 0) {
      add("citation_metadata_unverified", sprintf(
        "%d of %d references could not be fully verified (%d flagged, %d unverified, %d with warnings); see the audit array. Verification confirms identifiers and metadata only and does not check that each paper supports the specific claim.",
        bad, length(stat), sum(stat == "flagged"), sum(stat == "unverified"), sum(stat == "verified_with_warnings")),
        unique(unlist(lapply(audit[stat != "verified"], function(a) unlist(a$cited_by)))))
    }
  }

  # scope
  scope <- c(
    sprintf("Only %s were queried; no full-text literature, functional assays, or other databases.", paste(SOURCES[names(mods)], collapse = ", ")),
    "Coordinates are GRCh38; transcript-level annotation uses the single MANE/canonical transcript selected by Ensembl VEP.",
    "GTR evidence is gene-level and says nothing variant-specific.",
    "The score measures the weight and consistency of evidence, not whether the variant is pathogenic or benign."
  )
  if (!isTRUE(var$is_snv)) scope <- c(scope, "Variant is not a simple SNV: JASPAR scanning skipped and gnomAD queried by rsID.")
  if (is.na(var$protein_pos)) scope <- c(scope, "Variant has no protein position: AlphaFold/AlphaMissense not applicable.")
  if (isTRUE(mods$gwas$flags$truncated)) scope <- c(scope, "GWAS associations were truncated to the top traits and/or study look-ups were capped.")
  if (isTRUE(mods$clinvar$flags$truncated)) scope <- c(scope, sprintf("Only the first %d ClinVar records were evaluated.", CONFIG$clinvar_max_records))
  if (isTRUE((mods$clinvar$flags$scv_failed %||% 0) > 0)) scope <- c(scope, sprintf(
    "Submission-level ClinVar data were unavailable for %d record(s)%s, so their conflicts were scored from the aggregate classification only.",
    as.integer(mods$clinvar$flags$scv_failed),
    if (!requireNamespace("xml2", quietly = TRUE)) " (install the xml2 package to enable this)" else ""))
  if (length(var$alts) > 1) scope <- c(scope, sprintf("Multi-allelic site (%s); alleles are evaluated separately where sources allow.", var$allele_string))
  add("scope_restrictions", paste(scope, collapse = " "), srcs_of(names(mods)))

  # source availability
  failed <- names(st)[st == "error"]
  if (length(failed)) {
    add("source_unavailable", sprintf("Query failed for: %s. Missing sources may lower the score.",
                                      paste(vapply(failed, function(k) sprintf("%s (%s)", SOURCES[[k]], mods[[k]]$error %||% "error"), ""), collapse = "; ")),
        srcs_of(failed))
  }

  # data integrity
  di <- character()
  if (var$n_grch38_mappings > 1) di <- c(di, sprintf("rsID maps to %d GRCh38 locations; the first was used.", var$n_grch38_mappings))
  if (!identical(var$rsid, var$rsid_current)) di <- c(di, sprintf("Query %s was merged into %s; the current ID was used.", var$rsid, var$rsid_current))
  if (isTRUE(mods$alphafold$flags$isoform_mismatch)) di <- c(di, "Transcript residue differs from the UniProt sequence (possible isoform mismatch).")
  if (isTRUE(mods$jaspar$flags$ref_mismatch)) di <- c(di, "Genome reference base differs from the dbSNP reference allele.")
  if (isTRUE((mods$clinvar$flags$xref_unconfirmed %||% 0) > 0)) di <- c(di, "Some ClinVar records lacked a dbSNP cross-reference to confirm identity.")
  if (length(di)) add("data_integrity", paste(di, collapse = " "), srcs_of(c("alphafold", "jaspar", "clinvar")))

  L
}

# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------
direction_label <- function(d) {
  switch(as.character(d), "1" = "supports_functional_or_clinical_effect",
         "-1" = "argues_against_effect", "neutral_or_not_applicable")
}

rank_snp_evidence <- function(rsid, sources = names(SOURCES), out = NULL, verbose = TRUE) {
  rsid <- tolower(trimws(rsid))
  if (!grepl("^rs[0-9]+$", rsid)) stop("Query must be a dbSNP rsID such as 'rs429358'.", call. = FALSE)
  sources <- unique(match.arg(tolower(trimws(sources)), names(SOURCES), several.ok = TRUE))
  say <- function(...) if (verbose) message(...)
  retrieved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  say("Resolving ", rsid, " via Ensembl ...")
  var <- resolve_variant(rsid)

  runners <- list(clinvar = query_clinvar, gnomad = query_gnomad, gwas = query_gwas,
                  alphafold = query_alphafold, jaspar = query_jaspar, gtr = query_gtr)
  mods <- list()
  for (s in sources) {
    say("Querying ", SOURCES[[s]], " ...")
    mods[[s]] <- tryCatch(runners[[s]](var), error = function(e)
      module_result(SOURCES[[s]], "error", error = conditionMessage(e)))
  }

  items <- unname(do.call(c, lapply(mods, `[[`, "items")))
  xs <- cross_source_checks(mods)
  items <- c(items, xs$items)
  refs <- unname(do.call(c, lapply(mods, `[[`, "refs")))

  conflicts <- c(unname(do.call(c, lapply(mods, function(m) m$raw$conflicts %||% list()))), xs$conflicts)
  for (i in seq_along(conflicts)) conflicts[[i]]$id <- sprintf("C%d", i)
  ckeys <- vapply(conflicts, function(cf) paste(cf$source, cf$group, sep = "|"), "")
  for (i in seq_along(items)) {
    items[[i]]$id <- sprintf("E%d", i)
    k <- items[[i]]$conflict_key
    items[[i]]$conflict_id <- if (!is.na(k) && k %in% ckeys) conflicts[[match(k, ckeys)]]$id else NA_character_
  }

  say("Scoring and auditing ...")
  conf <- compute_confidence(items, conflicts)
  for (i in seq_along(conflicts)) {
    cid <- conflicts[[i]]$id
    conflicts[[i]]$booster_points <- round(unname(conf$booster_per_conflict[cid]), 2)
    conflicts[[i]]$affected_evidence_ids <- as.list(vapply(
      Filter(function(it) identical(it$conflict_id, cid), items), `[[`, "", "id"))
    conflicts[[i]]$effect_on_score <- if (identical(conflicts[[i]]$status, "resolved_correction"))
      "Superseded assertions carry no weight; the corrected position is scored as current evidence and a review booster is added."
    else "Recorded as mixed evidence: lowers the consistency factor and may trigger the Emerging gate."
  }
  conf$booster_per_conflict <- NULL
  audit <- build_audit(items, refs)
  limitations <- build_limitations(var, mods, items, audit, retrieved_at, conflicts)

  strength_rank <- match(vapply(items, `[[`, "", "strength"), c("direct", "strong", "indirect", "mixed", "emerging"))
  ord <- if (length(items)) order(strength_rank, -conf$item_points) else integer()
  evidence <- lapply(ord, function(i) {
    it <- items[[i]]
    list(id = it$id, source = it$source, strength = it$strength,
         strength_definition = STRENGTH_DEFS[[it$strength]], claim = it$claim,
         direction = direction_label(it$direction), replicated = it$replicated,
         quality_multiplier = round(it$quality, 2), weighted_points = round(conf$item_points[i], 2),
         conflict_context = it$conflict_context, conflict_id = it$conflict_id,
         independent_study_ids = as.list(it$study_ids), details = it$details,
         references = as.list(it$refs))
  })
  conf$item_points <- NULL

  var_out <- var
  var_out$alts <- as.list(var$alts); var_out$synonyms <- as.list(head(var$synonyms, 20))
  var_out$consequence_terms <- as.list(var$consequence_terms)

  result <- list(
    tool = list(name = "snp_evidence_ranker", version = TOOL_VERSION),
    query = list(rsid = rsid, sources = as.list(unname(SOURCES[sources]))),
    retrieved_at = retrieved_at,
    variant = var_out,
    source_status = lapply(mods, function(m) list(source = m$source, status = m$status,
                                                  error = m$error, n_items = length(m$items),
                                                  notes = as.list(m$notes))),
    confidence = conf,
    evidence_summary = evidence,
    conflicts = conflicts,
    limitations = limitations,
    audit = audit
  )
  if (!is.null(out)) {
    jsonlite::write_json(result, out, auto_unbox = TRUE, pretty = TRUE,
                         null = "null", na = "null", digits = NA)
    say("Wrote ", normalizePath(out, mustWork = FALSE))
  }
  invisible(result)
}

# -----------------------------------------------------------------------------
# Console report
# -----------------------------------------------------------------------------
print_report <- function(res) {
  w <- max(60, min(getOption("width", 100), 110))
  wrap <- function(txt, indent = 0, exdent = indent)
    cat(strwrap(txt, width = w - 2, indent = indent, exdent = exdent), sep = "\n")
  v <- res$variant
  cat("\n", strrep("=", w - 2), "\n", sep = "")
  cat(sprintf("SNP EVIDENCE REPORT: %s   (retrieved %s)\n", res$query$rsid, res$retrieved_at))
  cat(strrep("=", w - 2), "\n", sep = "")
  wrap(sprintf("Variant: chr%s:%s %s (GRCh38) | gene %s | %s%s",
               v$chr, format(v$start, scientific = FALSE), v$allele_string, v$gene %||% "NA",
               v$most_severe_consequence %||% "NA",
               if (length(v$missense)) paste0(" | p.", paste(vapply(v$missense, `[[`, "", "protein_change"), collapse = ", p.")) else ""))
  cat("\nSources: ")
  cat(paste(vapply(res$source_status, function(s) sprintf("%s [%s]", s$source, s$status), ""), collapse = ", "), "\n")

  cf <- res$confidence
  cat(sprintf("\nCONFIDENCE: %d / 100  ->  %s\n", cf$score, toupper(cf$band)))
  cat(sprintf("  baseline %.1f + correction booster %.1f (%d resolved, %d unresolved conflicts)\n",
              cf$baseline_score, cf$correction_booster,
              as.integer(cf$n_resolved_corrections), as.integer(cf$n_unresolved_conflicts)))
  wrap(cf$band_definition, 2)
  wrap(cf$rationale, 2)

  cat("\nEVIDENCE SUMMARY\n")
  if (!length(res$evidence_summary)) cat("  (no evidence items)\n")
  for (e in res$evidence_summary) {
    tag <- switch(e$conflict_context, corrected_by_review = " {CORRECTED}",
                  unresolved_conflict = " {UNRESOLVED}", "")
    wrap(sprintf("[%s] %-8s %s%s (%.1f pts): %s", e$id, toupper(e$strength), e$source, tag,
                 e$weighted_points, e$claim), 2, 6)
  }

  cat("\nCONFLICTS\n")
  if (!length(res$conflicts)) cat("  (none detected: all evidence scored at baseline)\n")
  for (cf in res$conflicts) {
    if (identical(cf$status, "resolved_correction")) {
      wrap(sprintf("%s RESOLVED (%s review, booster +%.1f) %s / %s: current position '%s' supersedes '%s'. Basis: %s.",
                   cf$id, cf$review_signal, cf$booster_points, cf$source, cf$group, cf$current_position,
                   paste(unlist(cf$superseded_positions), collapse = "/"),
                   paste(unlist(cf$signal_basis), collapse = "; ")), 2, 6)
    } else {
      wrap(sprintf("%s UNRESOLVED %s / %s: %s. %s", cf$id, cf$source, cf$group,
                   paste(unlist(cf$positions), collapse = " vs "), cf$reason %||% ""), 2, 6)
    }
  }

  cat("\nLIMITATIONS\n")
  for (l in res$limitations) wrap(sprintf("%s [%s] %s", l$id, l$code, l$description), 2, 6)

  cat("\nCITATION AUDIT\n")
  if (!length(res$audit)) {
    cat("  (no references)\n")
  } else {
    st <- vapply(res$audit, `[[`, "", "integrity_status")
    tab <- table(factor(st, levels = c("verified", "verified_with_warnings", "flagged", "unverified")))
    cat(sprintf("  %d references: %s\n", length(st), paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")))
    for (a in res$audit) {
      if (a$integrity_status %in% c("flagged", "unverified")) {
        wrap(sprintf("! %s (%s): %s", a$ref_id, a$integrity_status, paste(unlist(a$flags), collapse = ", ")), 4, 6)
      }
    }
  }
  cat(strrep("=", w - 2), "\n\n", sep = "")
  invisible(res)
}

# -----------------------------------------------------------------------------
# Command-line entry point
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste(
    "Usage: Rscript snp_evidence_ranker.R <rsID> [--out=FILE.json] [--sources=a,b,...] [--quiet]",
    sprintf("  sources: %s (default: all)", paste(names(SOURCES), collapse = ",")),
    "  env:     NCBI_API_KEY, NCBI_EMAIL (optional; faster NCBI queries)", sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) { cat(usage, "\n"); quit(status = 0) }
  opt <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
  }
  rsid_arg <- args[!startsWith(args, "--")][1]
  if (is.na(rsid_arg)) { cat(usage, "\n"); quit(status = 1) }
  out_file <- opt("out", paste0(tolower(rsid_arg), "_evidence.json"))
  src_arg <- strsplit(opt("sources", paste(names(SOURCES), collapse = ",")), ",")[[1]]
  res <- tryCatch(
    rank_snp_evidence(rsid_arg, sources = src_arg, out = out_file, verbose = !("--quiet" %in% args)),
    error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 1) })
  print_report(res)
}
