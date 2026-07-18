#!/usr/bin/env Rscript
# 02_process_sources.R
# Focused cleaning for Ensembl gene ID, gene symbol, and UniProt ID only.
# UniProt accession column is taken from Ensembl/HGNC cross-references
# (standalone UniProt download may be unavailable).

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}
source(file.path("scripts", "00_utils.R"))
ensure_project_dirs()

blank_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x)) | x == "-" | x == "NA"] <- NA_character_
  x
}

coalesce_chr <- function(...) {
  args <- list(...)
  n <- max(vapply(args, length, integer(1)))
  out <- rep(NA_character_, n)
  for (v in args) {
    v <- as.character(v)
    if (length(v) == 1L && n > 1L) v <- rep(v, n)
    pick <- is.na(out) | !nzchar(out)
    cand <- !is.na(v) & nzchar(v)
    out[pick & cand] <- v[pick & cand]
  }
  out
}

message("=== Step 02: Processing Ensembl / gene_symbol / UniProt only ===")

process_ensembl <- function() {
  message("\n--- Processing Ensembl BioMart (focused columns) ---")
  raw <- read_csv_safe(here_path("data", "raw", "ensembl_biomart_raw.csv"))
  message("Raw rows: ", nrow(raw))

  gene_symbol <- coalesce_chr(raw$external_gene_name, raw$hgnc_symbol)
  swiss <- blank_to_na(raw$uniprotswissprot)
  trembl <- blank_to_na(raw$uniprotsptrembl)
  uniprot <- ifelse(
    !is.na(swiss) & !is.na(trembl), paste(swiss, trembl, sep = ";"),
    ifelse(!is.na(swiss), swiss, trembl)
  )

  clean <- data.frame(
    ensembl_gene_id = blank_to_na(raw$ensembl_gene_id),
    gene_symbol     = blank_to_na(gene_symbol),
    uniprot_id      = blank_to_na(uniprot),
    source_database = "ensembl_biomart",
    mapping_source  = "ensembl_biomart",
    stringsAsFactors = FALSE
  )
  clean <- expand_all_multivalue_columns(clean, "uniprot_id")
  clean <- dedupe_rows_report(clean, "ensembl_biomart_clean")
  out <- here_path("data", "processed", "ensembl_biomart_clean.csv")
  write_csv_safe(clean, out)
  list(n_raw = nrow(raw), n_clean = nrow(clean), file = out)
}

process_hgnc <- function() {
  message("\n--- Processing HGNC (focused columns) ---")
  raw <- read_csv_safe(here_path("data", "raw", "hgnc_complete_raw.csv"))
  message("Raw rows: ", nrow(raw))

  clean <- data.frame(
    ensembl_gene_id = blank_to_na(raw$ensembl_gene_id),
    gene_symbol     = blank_to_na(raw$symbol),
    uniprot_id      = blank_to_na(raw$uniprot_ids),
    source_database = "hgnc",
    mapping_source  = "hgnc",
    stringsAsFactors = FALSE
  )
  clean <- expand_all_multivalue_columns(clean, c("ensembl_gene_id", "uniprot_id"))
  clean <- dedupe_rows_report(clean, "hgnc_complete_clean")
  out <- here_path("data", "processed", "hgnc_complete_clean.csv")
  write_csv_safe(clean, out)
  list(n_raw = nrow(raw), n_clean = nrow(clean), file = out)
}

stats <- list(
  ensembl = process_ensembl(),
  hgnc = process_hgnc()
)

# Optional UniProt raw if present; ignored when absent.
uniprot_raw <- here_path("data", "raw", "uniprot_human_raw.csv")
if (file.exists(uniprot_raw)) {
  message("\n--- UniProt raw present; processing focused columns ---")
  raw <- read_csv_safe(uniprot_raw)
  acc <- if ("uniprot_accession" %in% names(raw)) raw$uniprot_accession else raw$Entry
  gene <- if ("gene_names" %in% names(raw)) raw$gene_names else raw[["Gene Names"]]
  ens <- if ("xref_ensembl" %in% names(raw)) raw$xref_ensembl else raw$Ensembl
  primary <- vapply(as.character(gene), function(x) {
    p <- split_multivalue_string(x)
    if (length(p)) p[1] else NA_character_
  }, character(1))
  ens_gene <- vapply(as.character(ens), function(x) {
    vals <- split_multivalue_string(x)
    vals <- sub("\\s*\\[.*\\]$", "", vals)
    hit <- vals[grepl("^ENSG", vals)]
    if (length(hit)) paste(unique(hit), collapse = ";") else NA_character_
  }, character(1))
  clean <- data.frame(
    ensembl_gene_id = blank_to_na(ens_gene),
    gene_symbol     = blank_to_na(primary),
    uniprot_id      = blank_to_na(acc),
    source_database = "uniprot",
    mapping_source  = "uniprot",
    stringsAsFactors = FALSE
  )
  clean <- expand_all_multivalue_columns(clean, "ensembl_gene_id")
  clean <- dedupe_rows_report(clean, "uniprot_human_clean")
  out <- here_path("data", "processed", "uniprot_human_clean.csv")
  write_csv_safe(clean, out)
  stats$uniprot <- list(n_raw = nrow(raw), n_clean = nrow(clean), file = out)
} else {
  message("\nNote: UniProt raw download not found; using UniProt IDs from Ensembl + HGNC only.")
}

meta <- load_download_metadata()
meta$processing <- list(
  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  focus = c("ensembl_gene_id", "gene_symbol", "uniprot_id"),
  stats = stats
)
save_download_metadata(meta)

message("\n=== Processing complete ===")
