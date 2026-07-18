#!/usr/bin/env Rscript
# run_all.R — Focused pipeline (ensembl_gene_id / gene_symbol / uniprot_id).
# Assumes downloads already exist under data/raw/.

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}

steps <- c(
  "02_process_sources.R",
  "03_create_pairwise_mappings.R",
  "04_generate_unmapped_lists.R",
  "05_generate_report.R"
)

for (step in steps) {
  message("\n", strrep("=", 60))
  message("Running ", step)
  message(strrep("=", 60))
  status <- system2("Rscript", args = file.path("scripts", step))
  if (!identical(as.integer(status), 0L)) {
    stop("Pipeline failed at ", step, " (exit status ", status, ")")
  }
}

message("\nPipeline finished successfully.")
