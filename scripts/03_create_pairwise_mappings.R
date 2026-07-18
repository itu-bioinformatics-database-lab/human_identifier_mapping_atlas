#!/usr/bin/env Rscript
# 03_create_pairwise_mappings.R
# Pairwise mappings among ensembl_gene_id, gene_symbol, uniprot_id only.
# Outputs are JSON.

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}
source(file.path("scripts", "00_utils.R"))
ensure_project_dirs()

message("=== Step 03: Creating Ensembl / symbol / UniProt JSON mappings ===")

ensembl <- read_csv_safe(here_path("data", "processed", "ensembl_biomart_clean.csv"))
hgnc    <- read_csv_safe(here_path("data", "processed", "hgnc_complete_clean.csv"))
uniprot_path <- here_path("data", "processed", "uniprot_human_clean.csv")
uniprot <- if (file.exists(uniprot_path)) read_csv_safe(uniprot_path) else NULL

align_edge_cols <- function(a, b) {
  cols <- union(names(a), names(b))
  for (c in setdiff(cols, names(a))) a[[c]] <- NA
  for (c in setdiff(cols, names(b))) b[[c]] <- NA
  list(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}

PAIR_SPECS <- list(
  list(ensembl, "ensembl_gene_id", "gene_symbol", "ensembl_biomart", "ensembl_biomart"),
  list(ensembl, "ensembl_gene_id", "uniprot_id", "ensembl_biomart", "ensembl_biomart"),
  list(ensembl, "gene_symbol", "ensembl_gene_id", "ensembl_biomart", "ensembl_biomart"),
  list(ensembl, "gene_symbol", "uniprot_id", "ensembl_biomart", "ensembl_biomart"),
  list(ensembl, "uniprot_id", "ensembl_gene_id", "ensembl_biomart", "ensembl_biomart"),
  list(ensembl, "uniprot_id", "gene_symbol", "ensembl_biomart", "ensembl_biomart"),
  list(hgnc, "ensembl_gene_id", "gene_symbol", "hgnc", "hgnc"),
  list(hgnc, "ensembl_gene_id", "uniprot_id", "hgnc", "hgnc"),
  list(hgnc, "gene_symbol", "ensembl_gene_id", "hgnc", "hgnc"),
  list(hgnc, "gene_symbol", "uniprot_id", "hgnc", "hgnc"),
  list(hgnc, "uniprot_id", "ensembl_gene_id", "hgnc", "hgnc"),
  list(hgnc, "uniprot_id", "gene_symbol", "hgnc", "hgnc")
)

if (!is.null(uniprot)) {
  PAIR_SPECS <- c(
    PAIR_SPECS,
    list(
      list(uniprot, "ensembl_gene_id", "gene_symbol", "uniprot", "uniprot"),
      list(uniprot, "ensembl_gene_id", "uniprot_id", "uniprot", "uniprot"),
      list(uniprot, "gene_symbol", "ensembl_gene_id", "uniprot", "uniprot"),
      list(uniprot, "gene_symbol", "uniprot_id", "uniprot", "uniprot"),
      list(uniprot, "uniprot_id", "ensembl_gene_id", "uniprot", "uniprot"),
      list(uniprot, "uniprot_id", "gene_symbol", "uniprot", "uniprot")
    )
  )
}

NAMED_OUTPUTS <- c(
  "ensembl_gene_id.gene_symbol=ensembl_gene_to_gene_symbol.json",
  "ensembl_gene_id.uniprot_id=ensembl_gene_to_uniprot.json",
  "gene_symbol.ensembl_gene_id=gene_symbol_to_ensembl_gene.json",
  "gene_symbol.uniprot_id=gene_symbol_to_uniprot.json",
  "uniprot_id.ensembl_gene_id=uniprot_to_ensembl_gene.json",
  "uniprot_id.gene_symbol=uniprot_to_gene_symbol.json"
)

all_edges <- list()
edge_store <- list()

for (spec in PAIR_SPECS) {
  df <- spec[[1]]
  sc <- spec[[2]]
  tc <- spec[[3]]
  ms <- spec[[4]]
  sd <- spec[[5]]
  edges <- extract_pairwise_edges(
    df, sc, tc, ms, sd,
    context_cols = character(0)
  )
  if (nrow(edges) == 0) next
  all_edges[[length(all_edges) + 1]] <- edges
  key <- paste(sc, tc, sep = ".")
  if (is.null(edge_store[[key]])) {
    edge_store[[key]] <- edges
  } else {
    aligned <- align_edge_cols(edge_store[[key]], edges)
    edge_store[[key]] <- dedupe_rows_report(
      rbind(aligned[[1]], aligned[[2]]),
      key
    )
  }
}

for (entry in NAMED_OUTPUTS) {
  parts <- strsplit(entry, "=", fixed = TRUE)[[1]]
  key <- parts[1]
  fname <- parts[2]
  edges <- edge_store[[key]]
  if (is.null(edges)) {
    edges <- data.frame(
      source_id = character(),
      target_id = character(),
      source_id_type = character(),
      target_id_type = character(),
      mapping_source = character(),
      source_database = character(),
      stringsAsFactors = FALSE
    )
  }
  save_pairwise_table(edges, fname)
}

# Master long JSON
master_long <- combine_edges(all_edges)
master_long$confidence_note <- ""
keep_cols <- c(
  "source_id", "source_id_type", "target_id", "target_id_type",
  "mapping_source", "source_database", "confidence_note"
)
master_long <- master_long[, keep_cols, drop = FALSE]
write_json_safe(
  master_long,
  here_path("results", "mapping_tables", "master_identifier_mappings_long.json"),
  pretty = FALSE
)

# Master wide JSON: focused annotation rows
align_wide <- function(df) {
  std <- c("ensembl_gene_id", "gene_symbol", "uniprot_id", "source_database")
  out <- as.data.frame(
    matrix(NA_character_, nrow = nrow(df), ncol = length(std)),
    stringsAsFactors = FALSE
  )
  names(out) <- std
  for (col in intersect(std, names(df))) out[[col]] <- as.character(df[[col]])
  out
}
wide_parts <- list(align_wide(ensembl), align_wide(hgnc))
if (!is.null(uniprot)) wide_parts[[length(wide_parts) + 1]] <- align_wide(uniprot)
master_wide <- do.call(rbind, wide_parts)
rownames(master_wide) <- NULL
master_wide <- dedupe_rows_report(master_wide, "master_wide")
write_json_safe(
  master_wide,
  here_path("results", "mapping_tables", "master_identifier_annotations_wide.json"),
  pretty = FALSE
)

# Pairwise stats
stats_rows <- list()
for (key in names(edge_store)) {
  parts <- strsplit(key, ".", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    stats_rows[[length(stats_rows) + 1]] <- compute_pairwise_stats(
      edge_store[[key]], parts[1], parts[2]
    )
  }
}
pairwise_stats <- do.call(rbind, stats_rows)
rownames(pairwise_stats) <- NULL
write_csv_safe(pairwise_stats, here_path("results", "stats", "pairwise_mapping_stats.csv"))
write_json_safe(pairwise_stats, here_path("results", "stats", "pairwise_mapping_stats.json"))

count_unique <- function(df, col) {
  if (is.null(df) || !col %in% names(df)) return(0L)
  length(unique(df[[col]][!is.na(df[[col]]) & df[[col]] != ""]))
}

id_counts <- data.frame(
  id_type = c("ensembl_gene_id", "gene_symbol", "uniprot_id"),
  count = c(
    max(count_unique(ensembl, "ensembl_gene_id"), count_unique(hgnc, "ensembl_gene_id"),
        if (!is.null(uniprot)) count_unique(uniprot, "ensembl_gene_id") else 0L),
    max(count_unique(ensembl, "gene_symbol"), count_unique(hgnc, "gene_symbol"),
        if (!is.null(uniprot)) count_unique(uniprot, "gene_symbol") else 0L),
    max(count_unique(ensembl, "uniprot_id"), count_unique(hgnc, "uniprot_id"),
        if (!is.null(uniprot)) count_unique(uniprot, "uniprot_id") else 0L)
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(id_counts, here_path("results", "stats", "id_type_counts.csv"))
write_json_safe(id_counts, here_path("results", "stats", "id_type_counts.json"))

message("\n=== Pairwise mapping creation complete ===")
