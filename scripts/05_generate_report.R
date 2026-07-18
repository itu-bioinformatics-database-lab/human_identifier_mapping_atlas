#!/usr/bin/env Rscript
# 05_generate_report.R
# Focused report for ensembl_gene_id / gene_symbol / uniprot_id mappings.

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}
source(file.path("scripts", "00_utils.R"))
ensure_project_dirs()

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

message("=== Step 05: Generating focused markdown report ===")

meta <- load_download_metadata()
download_date <- meta$download_date %||% "unknown"
proc_date <- meta$processing$date %||% "unknown"

pairwise_stats <- read_csv_safe(here_path("results", "stats", "pairwise_mapping_stats.csv"))
unmapped_summary <- read_csv_safe(here_path("results", "stats", "unmapped_summary.csv"))
id_counts <- if (file.exists(here_path("results", "stats", "id_type_counts.csv"))) {
  read_csv_safe(here_path("results", "stats", "id_type_counts.csv"))
} else NULL

n_raw <- function(name) meta$processing$stats[[name]]$n_raw %||% NA_integer_
n_clean <- function(name) meta$processing$stats[[name]]$n_clean %||% NA_integer_
fmt_int <- function(x) if (is.na(x)) "N/A" else format(as.integer(x), big.mark = ",")

mapping_files <- sort(list.files(
  here_path("results", "mapping_tables"),
  pattern = "\\.json$",
  full.names = FALSE
))

lines <- c(
  "# Human Identifier Mapping Atlas — Report",
  "",
  paste("**Report generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 1. Project purpose",
  "",
  "Focused mapping atlas for three identifier types only:",
  "",
  "- Ensembl gene IDs (`ensembl_gene_id`)",
  "- Gene symbols (`gene_symbol`)",
  "- UniProt accessions (`uniprot_id`)",
  "",
  "All valid one-to-many / many-to-one / many-to-many mappings are preserved.",
  "All other ID types are ignored.",
  "",
  "## 2. Data sources used",
  "",
  "| Source | Status |",
  "|--------|--------|",
  "| Ensembl BioMart | Used (primary source) |",
  "| HGNC | Used (gene symbols + UniProt xrefs) |",
  paste0(
    "| UniProt proteome download | ",
    if (file.exists(here_path("data", "raw", "uniprot_human_raw.csv"))) {
      "Used"
    } else {
      "Not available — UniProt IDs taken from Ensembl/HGNC cross-references"
    },
    " |"
  ),
  "",
  "## 3. Download / processing dates",
  "",
  paste("- **Download metadata:**", download_date),
  paste("- **Processing completed:**", proc_date),
  "",
  "## 4. Source files retained",
  "",
  "- `data/raw/ensembl_biomart_raw.csv`",
  "- `data/raw/hgnc_complete_raw.csv`",
  if (file.exists(here_path("data", "raw", "uniprot_human_raw.csv"))) {
    "- `data/raw/uniprot_human_raw.csv`"
  } else {
    "- `data/raw/uniprot_human_raw.csv` *(missing)*"
  },
  "",
  "## 5. Cleaning steps",
  "",
  "1. Kept only `ensembl_gene_id`, `gene_symbol`, and `uniprot_id`.",
  "2. Split multi-value UniProt / Ensembl fields into long format.",
  "3. Removed only exact duplicate rows (counts reported in logs).",
  "4. Combined mappings from Ensembl and HGNC (and UniProt if present) with provenance.",
  "",
  "## 6. Raw / cleaned row counts",
  "",
  paste0("- **Ensembl raw:** ", fmt_int(n_raw("ensembl")),
         " → clean: ", fmt_int(n_clean("ensembl"))),
  paste0("- **HGNC raw:** ", fmt_int(n_raw("hgnc")),
         " → clean: ", fmt_int(n_clean("hgnc"))),
  if (!is.null(meta$processing$stats$uniprot)) {
    paste0("- **UniProt raw:** ", fmt_int(n_raw("uniprot")),
           " → clean: ", fmt_int(n_clean("uniprot")))
  } else {
    "- **UniProt raw:** not downloaded"
  },
  "",
  "## 7. Unique IDs by type",
  ""
)

if (!is.null(id_counts)) {
  lines <- c(lines, "| ID type | Count |", "|---------|------:|")
  for (i in seq_len(nrow(id_counts))) {
    lines <- c(lines, paste0("| ", id_counts$id_type[i], " | ", fmt_int(id_counts$count[i]), " |"))
  }
  lines <- c(lines, "")
}

lines <- c(lines,
  "## 8. Pairwise mapping counts",
  "",
  "| Source | Target | Unique sources | Unique targets | Mapping rows | 1:1 | 1:many | many:1 | many:many |",
  "|--------|--------|---------------:|---------------:|-------------:|----:|-------:|-------:|----------:|"
)
for (i in seq_len(nrow(pairwise_stats))) {
  r <- pairwise_stats[i, ]
  lines <- c(lines, paste0(
    "| ", r$source_id_type, " | ", r$target_id_type, " | ",
    fmt_int(r$n_unique_source_ids), " | ", fmt_int(r$n_unique_target_ids), " | ",
    fmt_int(r$n_mapping_rows), " | ", fmt_int(r$n_one_to_one), " | ",
    fmt_int(r$n_one_to_many), " | ", fmt_int(r$n_many_to_one), " | ",
    fmt_int(r$n_many_to_many), " |"
  ))
}

lines <- c(lines,
  "",
  "## 9. Mapping rates / unmapped counts",
  "",
  "| Source | Target | Total | Mapped | Unmapped | Mapping rate |",
  "|--------|--------|------:|-------:|---------:|-------------:|"
)
for (i in seq_len(nrow(unmapped_summary))) {
  r <- unmapped_summary[i, ]
  lines <- c(lines, paste0(
    "| ", r$source_id_type, " | ", r$target_id_type, " | ",
    fmt_int(r$total_source_ids), " | ", fmt_int(r$mapped_source_ids), " | ",
    fmt_int(r$unmapped_source_ids), " | ",
    if (is.na(r$mapping_rate)) "N/A" else sprintf("%.2f%%", r$mapping_rate * 100), " |"
  ))
}

lines <- c(lines,
  "",
  "## 10. Files generated",
  "",
  "### Mapping JSON (`results/mapping_tables/`)",
  ""
)
for (f in mapping_files) lines <- c(lines, paste0("- `", f, "`"))

lines <- c(lines,
  "",
  "### Stats",
  "",
  "- `results/stats/pairwise_mapping_stats.csv` / `.json`",
  "- `results/stats/unmapped_summary.csv` / `.json`",
  "- `results/stats/id_type_counts.csv` / `.json`",
  "",
  "### Unmapped lists (`data/unmapped/`)",
  "",
  paste0("- ", length(list.files(here_path("data", "unmapped"), pattern = "\\.txt$")), " text files"),
  "",
  "## 11. Known limitations",
  "",
  "- UniProt proteome download may be missing; UniProt IDs then come only from Ensembl/HGNC xrefs.",
  "- Not all Ensembl genes have UniProt IDs (especially non-coding genes).",
  "- Gene symbols can collide across loci; one-to-many mappings are expected and preserved.",
  "- Absence of a mapping does not mean the ID is invalid.",
  "",
  "---",
  "",
  "*Generated by `scripts/05_generate_report.R`*"
)

report_path <- here_path("results", "reports", "identifier_mapping_report.md")
writeLines(lines, report_path)
message("Wrote report: ", report_path)
message("\n=== Report generation complete ===")
