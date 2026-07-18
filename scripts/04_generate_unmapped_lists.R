#!/usr/bin/env Rscript
# 04_generate_unmapped_lists.R
# Unmapped lists for ensembl_gene_id, gene_symbol, uniprot_id only.

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}
source(file.path("scripts", "00_utils.R"))
ensure_project_dirs()

message("=== Step 04: Generating focused unmapped ID lists ===")

# Prefer JSON master; fall back to rebuilding from pairwise JSON if needed.
master_path_json <- here_path("results", "mapping_tables", "master_identifier_mappings_long.json")
if (!file.exists(master_path_json)) {
  stop("Required file not found: ", master_path_json)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package jsonlite is required to read master mappings JSON.")
}
master_long <- jsonlite::fromJSON(master_path_json, simplifyDataFrame = TRUE)

ensembl <- read_csv_safe(here_path("data", "processed", "ensembl_biomart_clean.csv"))
hgnc    <- read_csv_safe(here_path("data", "processed", "hgnc_complete_clean.csv"))
uniprot_path <- here_path("data", "processed", "uniprot_human_clean.csv")
uniprot <- if (file.exists(uniprot_path)) read_csv_safe(uniprot_path) else NULL

collect_ids <- function(id_type) {
  ids <- character(0)
  for (df in list(ensembl, hgnc, uniprot)) {
    if (!is.null(df) && id_type %in% names(df)) {
      vals <- df[[id_type]]
      vals <- vals[!is.na(vals) & vals != ""]
      ids <- c(ids, vals)
    }
  }
  unique(ids)
}

find_unmapped <- function(source_type, target_type) {
  all_source <- collect_ids(source_type)
  if (length(all_source) == 0) return(character(0))
  mapped <- master_long$source_id[
    master_long$source_id_type == source_type &
      master_long$target_id_type == target_type
  ]
  setdiff(all_source, unique(mapped))
}

write_unmapped_file <- function(ids, filename) {
  path <- here_path("data", "unmapped", filename)
  writeLines(sort(unique(ids)), path)
  message("Wrote ", length(unique(ids)), " unmapped IDs to ", filename)
  invisible(path)
}

UNMAPPED_SPECS <- list(
  c("ensembl_gene_id", "gene_symbol", "ensembl_genes_not_mapped_to_gene_symbol.txt"),
  c("ensembl_gene_id", "uniprot_id", "ensembl_genes_not_mapped_to_uniprot.txt"),
  c("gene_symbol", "ensembl_gene_id", "gene_symbols_not_mapped_to_ensembl_gene.txt"),
  c("gene_symbol", "uniprot_id", "gene_symbols_not_mapped_to_uniprot.txt"),
  c("uniprot_id", "ensembl_gene_id", "uniprots_not_mapped_to_ensembl_gene.txt"),
  c("uniprot_id", "gene_symbol", "uniprots_not_mapped_to_gene_symbol.txt")
)

summary_rows <- list()
for (spec in UNMAPPED_SPECS) {
  src_type <- spec[1]
  tgt_type <- spec[2]
  fname <- spec[3]
  all_src <- collect_ids(src_type)
  unmapped <- find_unmapped(src_type, tgt_type)
  write_unmapped_file(unmapped, fname)
  total <- length(all_src)
  n_unmapped <- length(unmapped)
  n_mapped <- total - n_unmapped
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    source_id_type = src_type,
    target_id_type = tgt_type,
    total_source_ids = total,
    mapped_source_ids = n_mapped,
    unmapped_source_ids = n_unmapped,
    mapping_rate = if (total > 0) round(n_mapped / total, 4) else NA_real_,
    unmapped_rate = if (total > 0) round(n_unmapped / total, 4) else NA_real_,
    unmapped_file = fname,
    stringsAsFactors = FALSE
  )
}

unmapped_summary <- do.call(rbind, summary_rows)
rownames(unmapped_summary) <- NULL
write_csv_safe(unmapped_summary, here_path("results", "stats", "unmapped_summary.csv"))
write_json_safe(unmapped_summary, here_path("results", "stats", "unmapped_summary.json"))

message("\n=== Unmapped list generation complete ===")
