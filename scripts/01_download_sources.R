#!/usr/bin/env Rscript
# 01_download_sources.R
# Download raw identifier tables from Ensembl BioMart, HGNC, and UniProt.
#
# BioMart is queried via the XML HTTP API (base R) so the pipeline does not
# require a working biomaRt binary package. If biomaRt is available it is used
# as a faster path.

args_cmd <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_cmd[grep("^--file=", args_cmd)])
if (length(script_path) == 1 && nzchar(script_path)) {
  setwd(dirname(normalizePath(script_path)))
  setwd("..")
}
source(file.path("scripts", "00_utils.R"))
ensure_project_dirs()

message("=== Step 01: Downloading raw data sources ===")
download_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
meta <- list(download_date = download_date, sources = list())

#' Merge BioMart data frames on shared key columns.
merge_on_keys <- function(a, b, keys) {
  keys <- intersect(keys, intersect(names(a), names(b)))
  if (length(keys) == 0) stop("No shared keys for BioMart merge.")
  merge(a, b, by = keys, all = TRUE, sort = FALSE)
}

#' Query Ensembl BioMart via XML HTTP POST (no biomaRt dependency).
biomart_query_http <- function(attributes, host = "https://www.ensembl.org",
                               dataset = "hsapiens_gene_ensembl") {
  attr_xml <- paste0(
    '      <Attribute name = "', attributes, '" />',
    collapse = "\n"
  )
  xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<!DOCTYPE Query>\n',
    '<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" ',
    'uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >\n',
    '      <Dataset name = "', dataset, '" interface = "default" >\n',
    attr_xml, "\n",
    "      </Dataset>\n",
    "</Query>"
  )
  url <- paste0(host, "/biomart/martservice")
  tmp <- tempfile(fileext = ".tsv")
  xml_file <- tempfile(fileext = ".xml")
  writeLines(xml, xml_file)
  ok <- FALSE
  # Write XML query to a file and POST via curl --data-urlencode query@file
  # to avoid shell newline quoting issues.
  status <- system2(
    "curl",
    c(
      "-sS", "-L",
      "--max-time", "1800",
      "-o", tmp,
      "-X", "POST",
      "--data-urlencode", paste0("query@", xml_file),
      url
    ),
    stdout = TRUE, stderr = TRUE
  )
  unlink(xml_file)
  if (file.exists(tmp) && file.info(tmp)$size >= 10) {
    ok <- TRUE
  } else {
    stop(
      "BioMart HTTP query failed for host ", host,
      if (length(status)) paste0(": ", paste(status, collapse = " ")) else ""
    )
  }
  first <- readLines(tmp, n = 5, warn = FALSE)
  if (any(grepl("ERROR|Exception|Query ERROR", first))) {
    unlink(tmp)
    stop(paste(first, collapse = " | "))
  }
  df <- utils::read.delim(tmp, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  unlink(tmp)
  # BioMart TSV headers are attribute display names; rename to attribute tokens if lengths match.
  if (ncol(df) == length(attributes)) {
    names(df) <- attributes
  } else {
    # attempt soft remapping when header names differ
    message("  Note: BioMart returned ", ncol(df), " columns for ", length(attributes), " attributes")
  }
  df
}

#' Download Ensembl using biomaRt package if it loads; otherwise use HTTP.
download_ensembl <- function() {
  message("\n--- Downloading Ensembl BioMart (hsapiens_gene_ensembl) ---")

  key_attrs <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id")
  query_sets <- list(
    core = c(
      key_attrs,
      "external_gene_name",
      "hgnc_symbol",
      "hgnc_id",
      "gene_biotype",
      "description"
    ),
    uniprot_entrez = c(key_attrs, "uniprotswissprot", "uniprotsptrembl", "entrezgene_id"),
    refseq = c(key_attrs, "refseq_mrna", "refseq_ncrna")
  )

  hosts <- c("https://www.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")
  used_host <- NULL
  use_biomart_pkg <- FALSE

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    biomart_ok <- tryCatch({
      suppressPackageStartupMessages(library(biomaRt))
      TRUE
    }, error = function(e) {
      message("biomaRt package present but failed to load: ", conditionMessage(e))
      message("Falling back to BioMart HTTP XML API.")
      FALSE
    })
    if (biomart_ok) {
      for (host in hosts) {
        message("Trying biomaRt host: ", host)
        mart <- tryCatch(
          biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = host),
          error = function(e) {
            message("  Failed: ", conditionMessage(e)); NULL
          }
        )
        if (!is.null(mart)) {
          used_host <- host
          use_biomart_pkg <- TRUE
          parts <- list()
          for (nm in names(query_sets)) {
            message("Fetching BioMart attribute set via biomaRt: ", nm, "...")
            parts[[nm]] <- biomaRt::getBM(attributes = query_sets[[nm]], mart = mart)
            message("  ", nm, ": ", nrow(parts[[nm]]), " rows")
          }
          raw <- parts$core
          raw <- merge_on_keys(raw, parts$uniprot_entrez, key_attrs)
          raw <- merge_on_keys(raw, parts$refseq, key_attrs)
          out_path <- here_path("data", "raw", "ensembl_biomart_raw.csv")
          write_csv_safe(raw, out_path)
          return(list(
            file = out_path, n_rows = nrow(raw), host = used_host,
            method = "biomaRt", n_cols = ncol(raw)
          ))
        }
      }
    }
  }

  # HTTP fallback
  raw <- NULL
  for (host in hosts) {
    message("Trying BioMart HTTP host: ", host)
    parts <- list()
    ok_host <- TRUE
    for (nm in names(query_sets)) {
      message("Fetching BioMart attribute set via HTTP: ", nm, "...")
      part <- tryCatch(
        biomart_query_http(query_sets[[nm]], host = host),
        error = function(e) {
          message("  Failed (", nm, "): ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(part)) {
        ok_host <- FALSE
        break
      }
      message("  ", nm, ": ", nrow(part), " rows")
      parts[[nm]] <- part
    }
    if (ok_host) {
      used_host <- host
      raw <- parts$core
      raw <- merge_on_keys(raw, parts$uniprot_entrez, key_attrs)
      raw <- merge_on_keys(raw, parts$refseq, key_attrs)
      break
    }
  }
  if (is.null(raw)) {
    stop("Could not download Ensembl BioMart from any host.")
  }

  message("Merged BioMart table: ", nrow(raw), " rows, ", ncol(raw), " columns")
  out_path <- here_path("data", "raw", "ensembl_biomart_raw.csv")
  write_csv_safe(raw, out_path)
  list(
    file = out_path, n_rows = nrow(raw), host = used_host,
    method = "http_xml", n_cols = ncol(raw)
  )
}

# --- 2. HGNC complete gene set ---
download_hgnc <- function() {
  message("\n--- Downloading HGNC complete gene set ---")
  urls <- c(
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_aliases&col=gd_prev_sym&col=gd_pub_ensembl_id&col=gd_pub_eg_id&col=gd_pub_uniprot_ids&col=gd_pub_refseq_ids&col=gd_locus_type&col=gd_status&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
  )
  out_path <- here_path("data", "raw", "hgnc_complete_raw.csv")

  raw <- NULL
  used_url <- NULL
  for (hgnc_url in urls) {
    message("Trying HGNC URL: ", hgnc_url)
    tmp <- tempfile(fileext = ".txt")
    dl_ok <- tryCatch({
      status <- system2("curl", c("-sS", "-L", "--max-time", "300", "-o", tmp, hgnc_url),
                        stdout = TRUE, stderr = TRUE)
      file.exists(tmp) && file.info(tmp)$size > 1000
    }, error = function(e) {
      message("  curl failed: ", conditionMessage(e))
      FALSE
    })
    if (!dl_ok) {
      dl_ok <- tryCatch({
        utils::download.file(hgnc_url, tmp, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) {
        message("  Download failed: ", conditionMessage(e))
        FALSE
      })
    }
    if (!dl_ok || !file.exists(tmp) || file.info(tmp)$size < 1000) {
      unlink(tmp)
      next
    }
    raw <- tryCatch(
      utils::read.delim(tmp, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) {
        message("  Parse failed: ", conditionMessage(e))
        NULL
      }
    )
    unlink(tmp)
    if (!is.null(raw) && nrow(raw) > 0) {
      used_url <- hgnc_url
      break
    }
  }
  if (is.null(raw)) stop("HGNC download failed from all candidate URLs.")

  rename_map <- c(
    "HGNC ID" = "hgnc_id",
    "Approved symbol" = "symbol",
    "Approved name" = "name",
    "Alias symbols" = "alias_symbol",
    "Previous symbols" = "prev_symbol",
    "Ensembl ID(supplied by Ensembl)" = "ensembl_gene_id",
    "NCBI Gene ID(supplied by NCBI)" = "entrez_id",
    "UniProt ID(supplied by UniProt)" = "uniprot_ids",
    "RefSeq(supplied by NCBI)" = "refseq_accession",
    "Locus type" = "locus_type",
    "Status" = "status"
  )
  for (old in names(rename_map)) {
    if (old %in% names(raw) && !rename_map[[old]] %in% names(raw)) {
      names(raw)[names(raw) == old] <- rename_map[[old]]
    }
  }

  keep <- c(
    "hgnc_id", "symbol", "name", "alias_symbol", "prev_symbol",
    "ensembl_gene_id", "entrez_id", "uniprot_ids", "refseq_accession",
    "locus_type", "status"
  )
  keep <- intersect(keep, names(raw))
  if (!"hgnc_id" %in% keep || !"symbol" %in% keep) {
    stop("HGNC file missing required columns. Found: ", paste(names(raw), collapse = ", "))
  }
  raw <- raw[, keep, drop = FALSE]
  write_csv_safe(raw, out_path)
  list(file = out_path, n_rows = nrow(raw), url = used_url)
}

# --- 3. UniProt human proteome via streaming endpoint ---
download_uniprot <- function() {
  message("\n--- Downloading UniProt human proteome (organism_id:9606) ---")
  fields <- paste(c(
    "accession", "reviewed", "protein_name", "gene_names",
    "xref_ensembl", "xref_hgnc", "xref_geneid", "xref_refseq", "organism_id"
  ), collapse = ",")

  out_path <- here_path("data", "raw", "uniprot_human_raw.csv")
  tmp_tsv <- tempfile(fileext = ".tsv")
  stream_url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "format=tsv&",
    "query=", utils::URLencode("organism_id:9606", reserved = TRUE), "&",
    "fields=", utils::URLencode(fields, reserved = TRUE)
  )

  message("Streaming UniProt TSV...")
  dl_ok <- FALSE
  # Prefer curl (handles large streams more reliably than download.file)
  curl_status <- system2(
    "curl",
    c("-sS", "-L", "--max-time", "1800", "-o", tmp_tsv, stream_url),
    stdout = TRUE, stderr = TRUE
  )
  if (file.exists(tmp_tsv) && file.info(tmp_tsv)$size > 1000) {
    dl_ok <- TRUE
  } else {
    message("curl stream failed: ", paste(curl_status, collapse = " "))
    status <- tryCatch({
      utils::download.file(stream_url, tmp_tsv, mode = "wb", quiet = FALSE)
      0
    }, error = function(e) {
      message("download.file failed: ", conditionMessage(e))
      1
    })
    dl_ok <- identical(as.integer(status), 0L) && file.exists(tmp_tsv) && file.info(tmp_tsv)$size > 1000
  }
  if (!dl_ok) {
    unlink(tmp_tsv)
    stop("UniProt stream download failed.")
  }

  message("Parsing streamed UniProt file (", round(file.info(tmp_tsv)$size / 1e6, 1), " MB)...")
  raw <- utils::read.delim(
    tmp_tsv,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  unlink(tmp_tsv)

  rename_exact <- c(
    "Entry" = "uniprot_accession",
    "Reviewed" = "reviewed_status",
    "Protein names" = "protein_name",
    "Gene Names" = "gene_names",
    "Ensembl" = "xref_ensembl",
    "HGNC" = "xref_hgnc",
    "GeneID" = "xref_geneid",
    "RefSeq" = "xref_refseq",
    "Organism (ID)" = "organism_id"
  )
  for (old in names(rename_exact)) {
    if (old %in% names(raw)) names(raw)[names(raw) == old] <- rename_exact[[old]]
  }

  if (!"uniprot_accession" %in% names(raw)) {
    stop("UniProt download missing accession column. Columns: ", paste(names(raw), collapse = ", "))
  }
  if ("organism_id" %in% names(raw)) {
    raw <- raw[as.character(raw$organism_id) == "9606", , drop = FALSE]
  }

  write_csv_safe(raw, out_path)
  list(
    file = out_path,
    n_rows = nrow(raw),
    url = "https://rest.uniprot.org/uniprotkb/stream",
    query = "organism_id:9606"
  )
}

force <- identical(Sys.getenv("FORCE_REDOWNLOAD", ""), "1")

existing_ok <- function(path, min_bytes = 1000) {
  if (force || !file.exists(path)) return(FALSE)
  info <- file.info(path)
  !is.na(info$size) && info$size >= min_bytes
}

run_download <- function(name, fn, out_file) {
  if (existing_ok(out_file)) {
    message("Skipping ", name, " (existing file): ", out_file)
    n <- tryCatch({
      # approximate row count without full load
      length(readLines(out_file, warn = FALSE)) - 1L
    }, error = function(e) NA_integer_)
    meta$sources[[name]] <<- list(file = out_file, n_rows = n, skipped = TRUE)
    return(TRUE)
  }
  tryCatch({
    result <- fn()
    meta$sources[[name]] <<- result
    message("SUCCESS: ", name)
    TRUE
  }, error = function(e) {
    message("ERROR downloading ", name, ": ", conditionMessage(e))
    meta$sources[[name]] <<- list(error = conditionMessage(e))
    FALSE
  })
}

ok_ensembl <- run_download(
  "ensembl_biomart", download_ensembl,
  here_path("data", "raw", "ensembl_biomart_raw.csv")
)
ok_hgnc <- run_download(
  "hgnc", download_hgnc,
  here_path("data", "raw", "hgnc_complete_raw.csv")
)
ok_uniprot <- run_download(
  "uniprot", download_uniprot,
  here_path("data", "raw", "uniprot_human_raw.csv")
)

meta$download_date <- download_date
save_download_metadata(meta)

if (!all(c(ok_ensembl, ok_hgnc, ok_uniprot))) {
  message("\nSome downloads failed. Re-run this script when network is available.")
  message("Set FORCE_REDOWNLOAD=1 to re-download existing files.")
  quit(status = 1)
}

message("\n=== Download complete ===")
