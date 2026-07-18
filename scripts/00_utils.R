# Shared utilities for human_identifier_mapping_atlas
# Sourced by all pipeline scripts. Base-R only (no hard CRAN deps).

#' Resolve project root from the current working directory.
project_root <- local({
  .cached <- NULL
  function(reset = FALSE) {
    if (!reset && !is.null(.cached)) return(.cached)
    candidates <- c(
      getwd(),
      if (basename(getwd()) == "scripts") dirname(getwd()) else NA_character_
    )
    root <- NULL
    for (p in unique(na.omit(candidates))) {
      if (file.exists(file.path(p, "README.md")) && dir.exists(file.path(p, "scripts"))) {
        root <- normalizePath(p, winslash = "/", mustWork = TRUE)
        break
      }
    }
    if (is.null(root)) root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    .cached <<- root
    root
  }
})

#' Build a path relative to the project root.
here_path <- function(...) {
  file.path(project_root(), ...)
}

# Prefer user library when available (optional packages like jsonlite).
local({
  home <- path.expand("~")
  for (lib in c(
    file.path(home, "R", "x86_64-pc-linux-gnu-library", "4.6"),
    file.path(home, "R", "x86_64-pc-linux-gnu-library", "4.5"),
    file.path(project_root(), ".R_libs")
  )) {
    if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  }
})

# Paths use here_path() (project-root relative). No required CRAN deps.

#' Create project directories if missing.
ensure_project_dirs <- function() {
  dirs <- c(
    "data/raw",
    "data/processed",
    "data/unmapped",
    "results/mapping_tables",
    "results/stats",
    "results/reports"
  )
  for (d in dirs) {
    path <- here_path(d)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      message("Created directory: ", path)
    }
  }
  invisible(NULL)
}

#' Split a single string on common multi-value delimiters.
split_multivalue_string <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "" || x == "-") {
    return(character(0))
  }
  x <- trimws(as.character(x))
  parts <- unlist(strsplit(x, "[,|;\\s]+", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & parts != "-" & parts != "NA"]
  unique(parts)
}

#' Expand one column with multi-value delimiters into long format.
expand_multivalue_column <- function(df, col_name) {
  if (!col_name %in% names(df) || nrow(df) == 0) {
    return(df)
  }
  vals <- as.character(df[[col_name]])
  needs_split <- !is.na(vals) & grepl("[,|;[:space:]]", vals)
  if (!any(needs_split)) {
    return(df)
  }
  split_list <- lapply(vals, split_multivalue_string)
  split_list <- lapply(seq_along(split_list), function(i) {
    if (length(split_list[[i]]) == 0) {
      if (is.na(vals[i]) || !nzchar(vals[i])) return(NA_character_)
      return(vals[i])
    }
    split_list[[i]]
  })
  n_rep <- lengths(split_list)
  out <- df[rep(seq_len(nrow(df)), n_rep), , drop = FALSE]
  out[[col_name]] <- unlist(split_list, use.names = FALSE)
  rownames(out) <- NULL
  out
}

#' Sequentially expand all listed multi-value columns.
expand_all_multivalue_columns <- function(df, columns) {
  out <- df
  for (col in columns) {
    if (col %in% names(out)) {
      before <- nrow(out)
      out <- expand_multivalue_column(out, col)
      after <- nrow(out)
      if (after > before) {
        message("  Expanded ", col, ": ", before, " -> ", after, " rows")
      }
    }
  }
  out
}

#' Remove exact duplicate rows and report how many were removed.
dedupe_rows_report <- function(df, label = "table") {
  before <- nrow(df)
  out <- unique(df)
  after <- nrow(out)
  removed <- before - after
  if (removed > 0) {
    message(
      "Removed ", removed, " exact duplicate row(s) from ", label,
      " (", before, " -> ", after, ")"
    )
  }
  out
}

#' Write CSV with consistent settings.
write_csv_safe <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
  message("Wrote: ", path, " (", nrow(df), " rows)")
}

#' Read CSV safely.
read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Extract pairwise edges between two ID columns from a data frame.
extract_pairwise_edges <- function(
    df,
    source_col,
    target_col,
    mapping_source,
    source_database,
    context_cols = c("gene_symbol", "hgnc_id", "gene_biotype", "reviewed_status", "description")
) {
  if (!source_col %in% names(df) || !target_col %in% names(df)) {
    return(data.frame())
  }
  sub <- df[!is.na(df[[source_col]]) & df[[source_col]] != "" &
            !is.na(df[[target_col]]) & df[[target_col]] != "", , drop = FALSE]
  if (nrow(sub) == 0) {
    return(data.frame())
  }
  out <- data.frame(
    source_id         = sub[[source_col]],
    target_id         = sub[[target_col]],
    source_id_type    = source_col,
    target_id_type    = target_col,
    mapping_source    = mapping_source,
    source_database   = source_database,
    stringsAsFactors  = FALSE
  )
  for (cc in context_cols) {
    if (cc %in% names(sub)) {
      out[[cc]] <- sub[[cc]]
    }
  }
  dedupe_rows_report(out, paste0(source_col, "_to_", target_col))
}

#' Combine edges from multiple data frames.
combine_edges <- function(edge_list) {
  if (length(edge_list) == 0) {
    return(data.frame())
  }
  edge_list <- edge_list[vapply(edge_list, nrow, integer(1)) > 0]
  if (length(edge_list) == 0) {
    return(data.frame())
  }
  all_cols <- unique(unlist(lapply(edge_list, names)))
  edge_list <- lapply(edge_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (m in missing) df[[m]] <- NA
    df[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, edge_list)
  rownames(out) <- NULL
  dedupe_rows_report(out, "combined edges")
}

#' Compute pairwise mapping cardinality stats for one edge table.
compute_pairwise_stats <- function(edges, source_type, target_type) {
  if (nrow(edges) == 0) {
    return(data.frame(
      source_id_type = source_type,
      target_id_type = target_type,
      n_unique_source_ids = 0L,
      n_unique_target_ids = 0L,
      n_mapping_rows = 0L,
      n_one_to_one = 0L,
      n_one_to_many = 0L,
      n_many_to_one = 0L,
      n_many_to_many = 0L,
      mapping_rate = NA_real_,
      source_database_counts = "",
      stringsAsFactors = FALSE
    ))
  }
  # Unique source->target and target->source degree maps
  src_n_tgt <- as.integer(tapply(edges$target_id, edges$source_id, function(x) length(unique(x))))
  tgt_n_src <- as.integer(tapply(edges$source_id, edges$target_id, function(x) length(unique(x))))
  names(src_n_tgt) <- names(tapply(edges$target_id, edges$source_id, length))
  names(tgt_n_src) <- names(tapply(edges$source_id, edges$target_id, length))
  n_src <- length(src_n_tgt)
  n_tgt <- length(tgt_n_src)
  n_rows <- nrow(edges)

  sc <- src_n_tgt[as.character(edges$source_id)]
  tc <- tgt_n_src[as.character(edges$target_id)]
  n_one_to_one <- sum(sc == 1 & tc == 1, na.rm = TRUE)
  n_one_to_many <- sum(sc > 1 & tc == 1, na.rm = TRUE)
  n_many_to_one <- sum(sc == 1 & tc > 1, na.rm = TRUE)
  n_many_to_many <- sum(sc > 1 & tc > 1, na.rm = TRUE)

  db_counts <- if ("source_database" %in% names(edges)) {
    paste(names(table(edges$source_database)), table(edges$source_database), sep = ":", collapse = ";")
  } else {
    ""
  }
  data.frame(
    source_id_type = source_type,
    target_id_type = target_type,
    n_unique_source_ids = n_src,
    n_unique_target_ids = n_tgt,
    n_mapping_rows = n_rows,
    n_one_to_one = n_one_to_one,
    n_one_to_many = n_one_to_many,
    n_many_to_one = n_many_to_one,
    n_many_to_many = n_many_to_many,
    mapping_rate = round(n_rows / max(n_src, 1), 4),
    source_database_counts = db_counts,
    stringsAsFactors = FALSE
  )
}

#' Write JSON (uses jsonlite when available; otherwise base R JSON encoder).
write_json_safe <- function(x, path, pretty = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path, pretty = pretty, auto_unbox = TRUE, dataframe = "rows", na = "null")
  } else {
    # Minimal fallback for list-of-named-lists / data.frame rows
    to_json_fallback <- function(obj) {
      if (is.null(obj)) return("null")
      if (is.logical(obj) && length(obj) == 1) return(if (is.na(obj)) "null" else tolower(as.character(obj)))
      if (is.numeric(obj) && length(obj) == 1) return(if (is.na(obj)) "null" else as.character(obj))
      if (is.character(obj) && length(obj) == 1) {
        if (is.na(obj)) return("null")
        esc <- gsub("\\", "\\\\", obj, fixed = TRUE)
        esc <- gsub("\"", "\\\"", esc, fixed = TRUE)
        esc <- gsub("\n", "\\n", esc, fixed = TRUE)
        esc <- gsub("\r", "\\r", esc, fixed = TRUE)
        esc <- gsub("\t", "\\t", esc, fixed = TRUE)
        return(paste0("\"", esc, "\""))
      }
      if (is.data.frame(obj)) {
        rows <- lapply(seq_len(nrow(obj)), function(i) {
          vals <- lapply(names(obj), function(nm) {
            paste0("\"", nm, "\":", to_json_fallback(obj[[nm]][i]))
          })
          paste0("{", paste(vals, collapse = ","), "}")
        })
        return(paste0("[", paste(rows, collapse = ","), "]"))
      }
      if (is.list(obj)) {
        if (is.null(names(obj))) {
          return(paste0("[", paste(vapply(obj, to_json_fallback, character(1)), collapse = ","), "]"))
        }
        parts <- vapply(names(obj), function(nm) {
          paste0("\"", nm, "\":", to_json_fallback(obj[[nm]]))
        }, character(1))
        return(paste0("{", paste(parts, collapse = ","), "}"))
      }
      if (length(obj) > 1) {
        return(paste0("[", paste(vapply(obj, to_json_fallback, character(1)), collapse = ","), "]"))
      }
      "null"
    }
    txt <- to_json_fallback(x)
    writeLines(txt, path)
  }
  n <- if (is.data.frame(x)) nrow(x) else if (is.list(x)) length(x) else NA
  message("Wrote: ", path, if (!is.na(n)) paste0(" (", n, " records)") else "")
  invisible(path)
}

#' Save a named pairwise mapping file as JSON.
save_pairwise_table <- function(edges, filename, pretty = FALSE) {
  if (!grepl("\\.json$", filename)) {
    filename <- sub("\\.csv$", ".json", filename)
    if (!grepl("\\.json$", filename)) filename <- paste0(filename, ".json")
  }
  if (nrow(edges) == 0) {
    message("Warning: empty mapping table for ", filename, " — writing empty array")
  }
  path <- here_path("results", "mapping_tables", filename)
  write_json_safe(edges, path, pretty = pretty)
  invisible(path)
}

#' Get download metadata path.
metadata_path <- function() {
  here_path("data", "raw", "download_metadata.rds")
}

#' Save download metadata.
save_download_metadata <- function(meta) {
  saveRDS(meta, metadata_path())
}

#' Load download metadata.
load_download_metadata <- function() {
  if (file.exists(metadata_path())) {
    readRDS(metadata_path())
  } else {
    list()
  }
}
