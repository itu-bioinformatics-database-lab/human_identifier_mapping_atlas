# Data Sources & Download Guide

This document explains **exactly how the raw data for each database is downloaded**,
where it lands on disk, how it is processed, and which files are intentionally
**excluded from git** because of their size.

The atlas maps between these human identifier types:

- `ensembl_gene_id` (e.g. `ENSG00000141510`)
- `ensembl_transcript_id` (e.g. `ENST00000361390`) — from BioMart raw via `scripts/07_create_transcript_mappings.py`
- `gene_symbol` (HGNC approved symbol, e.g. `TP53`)
- `uniprot_id` (UniProt accession, e.g. `P04637`)

All one-to-many / many-to-one / many-to-many relationships are preserved.

Transcript pairwise tables (same JSON edge schema as gene/symbol/UniProt):

- `ensembl_transcript_to_ensembl_gene.json`
- `ensembl_transcript_to_uniprot.json`
- `ensembl_transcript_to_gene_symbol.json`

---

## 1. Overview of the download step

All raw downloads are driven by **`scripts/01_download_sources.R`**. Run it with:

```bash
cd ~/Desktop/human_identifier_mapping_atlas
Rscript scripts/01_download_sources.R
# force re-download of files that already exist:
FORCE_REDOWNLOAD=1 Rscript scripts/01_download_sources.R
```

The script:

1. downloads **three** sources (Ensembl BioMart, HGNC, UniProt),
2. writes each raw table into `data/raw/`,
3. records a provenance record (date, host/URL, row counts, method) into
   `data/raw/download_metadata.rds`,
4. skips any source whose output already exists (unless `FORCE_REDOWNLOAD=1`).

Only the system tool **`curl`** is required. No CRAN packages are mandatory —
`biomaRt` is used automatically *if* it loads, otherwise the script falls back to a
plain HTTP call.

| Source | Raw output file | Transport | Approx size |
|--------|-----------------|-----------|-------------|
| Ensembl BioMart | `data/raw/ensembl_biomart_raw.csv` | biomaRt pkg → HTTP XML fallback | 50–150 MB |
| HGNC | `data/raw/hgnc_complete_raw.csv` | `curl` → `download.file` fallback | 5–10 MB |
| UniProt | `data/raw/uniprot_human_raw.csv` | `curl` stream → `download.file` fallback | 100–300 MB |

---

## 2. Ensembl BioMart

**Function:** `download_ensembl()` in `scripts/01_download_sources.R`
**Dataset:** `hsapiens_gene_ensembl`

### How it downloads

Two paths, tried in order:

1. **biomaRt package (fast path)** — if `biomaRt` is installed *and* loads, the script
   calls `biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = ...)`
   and fetches attributes with `getBM()`.
2. **HTTP XML fallback** — if biomaRt is missing or fails to load (e.g. an R 4.6 runtime
   with packages compiled for R 4.5), the script builds a BioMart XML query and POSTs it
   with `curl` to `<host>/biomart/martservice` (`--data-urlencode query@<file>`, TSV
   formatter). This needs no Bioconductor install.

**Host failover** (tried in order):

1. `https://www.ensembl.org`
2. `https://useast.ensembl.org`
3. `https://asia.ensembl.org`

### Split queries (why three requests)

BioMart rejects too many external-reference attributes in one query
(*"Too many attributes selected for External References"*). To avoid this, attributes are
fetched in **three separate queries keyed on the Ensembl gene/transcript/peptide IDs** and
then merged on those keys:

| Query set | Attributes |
|-----------|------------|
| `core` | `ensembl_gene_id`, `ensembl_transcript_id`, `ensembl_peptide_id`, `external_gene_name`, `hgnc_symbol`, `hgnc_id`, `gene_biotype`, `description` |
| `uniprot_entrez` | key IDs + `uniprotswissprot`, `uniprotsptrembl`, `entrezgene_id` |
| `refseq` | key IDs + `refseq_mrna`, `refseq_ncrna` |

The merged table is written to **`data/raw/ensembl_biomart_raw.csv`**.

---

## 3. HGNC (HUGO Gene Nomenclature Committee)

**Function:** `download_hgnc()`
**Format:** TSV

### How it downloads

`curl` is tried first, then `utils::download.file()` as a fallback. Two URLs are attempted
in order:

1. **Primary (complete set):**
   `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
2. **Fallback (genenames.org custom download):** a `cgi-bin/download/custom?...` URL that
   requests HGNC ID, approved symbol/name, aliases, previous symbols, Ensembl ID, Entrez ID,
   UniProt IDs, RefSeq, locus type, and status.

### Column normalization

Human-readable headers are renamed to canonical tokens, then only the useful columns are
kept:

`hgnc_id`, `symbol`, `name`, `alias_symbol`, `prev_symbol`, `ensembl_gene_id`, `entrez_id`,
`uniprot_ids`, `refseq_accession`, `locus_type`, `status`.

If `hgnc_id` or `symbol` is missing the download is treated as invalid. The cleaned table is
written to **`data/raw/hgnc_complete_raw.csv`**.

---

## 4. UniProt (human proteome)

**Function:** `download_uniprot()`
**Endpoint:** `https://rest.uniprot.org/uniprotkb/stream` (TSV)
**Query:** `organism_id:9606` (reviewed Swiss-Prot **and** unreviewed TrEMBL)

### How it downloads

The script builds a streaming REST URL requesting these fields:

```
accession, reviewed, protein_name, gene_names,
xref_ensembl, xref_hgnc, xref_geneid, xref_refseq, organism_id
```

It streams with `curl` (`--max-time 1800`), and falls back to `utils::download.file()` if
curl fails. After download it:

- renames columns (`Entry → uniprot_accession`, `Gene Names → gene_names`,
  `Ensembl → xref_ensembl`, etc.),
- re-filters to `organism_id == 9606`,
- writes **`data/raw/uniprot_human_raw.csv`**.

> **If the UniProt download fails**, the pipeline still works: UniProt accessions are then
> taken from the Ensembl BioMart and HGNC cross-reference columns instead. `02_process_sources.R`
> treats `uniprot_human_raw.csv` as optional.

---

## 5. Artifact mappers (pre-downloaded UniProt / HGNC dumps)

Independent of the R downloader, you can drop **already-downloaded** source dumps under
`artifacts/mappers/` and build per-source mapping tables with a standalone Python script
(`scripts/06_process_artifact_mappers.py`, base Python only):

```text
artifacts/mappers/UniProt/HUMAN_9606_idmapping.dat.gz   # from UniProt FTP/idmapping
artifacts/mappers/HGNC/hgnc_complete_set.txt            # HGNC complete set TSV
```

```bash
python3 scripts/06_process_artifact_mappers.py
```

Where these dumps come from:

- **UniProt idmapping:** `HUMAN_9606_idmapping.dat.gz` — the human ID-mapping file from
  UniProt's downloads (`knowledgebase/idmapping/by_organism/`). Long format:
  `accession <TAB> id_type <TAB> value`. The script reads `Gene_Name` and `Ensembl`
  records, strips Ensembl version suffixes, and collapses isoforms (`P31946-1 → P31946`).
- **HGNC complete set:** same `hgnc_complete_set.txt` as the primary HGNC URL above.

For each source it writes, **into that source's own folder**:
`mappings/*.json`, `unmapped/*.txt`, `stats/*.csv|json`, and
`reports/identifier_mapping_report.md`. Raw download files are never deleted.

---

## 6. Processing pipeline (after download)

Run everything (assumes `data/raw/` is populated):

```bash
Rscript scripts/run_all.R
```

Which runs, in order:

| Step | Script | Purpose |
|------|--------|---------|
| 02 | `02_process_sources.R` | Keep only the 3 ID columns per source; split multi-value UniProt/Ensembl fields into long format; dedupe → `data/processed/*_clean.csv` |
| 03 | `03_create_pairwise_mappings.R` | Build the 6 directed pairwise JSON tables + master long/wide JSON; compute cardinality stats |
| 04 | `04_generate_unmapped_lists.R` | Per direction, list source IDs with no mapping → `data/unmapped/*.txt` + `unmapped_summary` |
| 05 | `05_generate_report.R` | Assemble `results/reports/identifier_mapping_report.md` |

Shared helpers (path resolution, multi-value splitting, dedupe, JSON writing, stats) live in
`scripts/00_utils.R`.

### Output files

| Path | Description |
|------|-------------|
| `results/mapping_tables/*.json` | 6 pairwise tables + `master_identifier_mappings_long.json` + `master_identifier_annotations_wide.json` |
| `results/stats/` | pairwise stats, ID-type counts, unmapped summary |
| `data/unmapped/*.txt` | one unmapped ID per line, per direction |
| `results/reports/identifier_mapping_report.md` | human-readable summary |

---

## 7. Large files & what git ignores

Raw downloads and generated mapping tables are **large and fully regenerable**, so they are
excluded from version control by `.gitignore`. The biggest offenders exceed GitHub's 100 MB
per-file limit and must never be committed:

| File | Size | Tracked? |
|------|------|----------|
| `artifacts/mappers/UniProt/mappings/master_identifier_mappings_long.json` | ~126 MB | ignored |
| `results/mapping_tables/master_identifier_mappings_long.json` | ~120 MB | ignored |
| `data/raw/ensembl_biomart_raw.csv` | ~117 MB | ignored |
| `artifacts/mappers/UniProt/HUMAN_9606_idmapping.dat.gz` | ~37 MB | ignored |
| `artifacts/mappers/HGNC/hgnc_complete_set.txt` | ~17 MB | ignored |
| all other `mappings/*.json`, `data/raw/*`, `data/processed/*` | MB-scale | ignored |

**What stays in the repo:** all `scripts/`, `README.md`, `requirements_notes.md`,
`DATA_SOURCES.md`, and the small generated `reports/*.md`.

To rebuild everything from scratch after cloning:

```bash
Rscript scripts/01_download_sources.R   # re-download raw sources
Rscript scripts/run_all.R               # process → mappings → unmapped → report
# optional, if you have the artifact dumps:
python3 scripts/06_process_artifact_mappers.py
```

---

## 8. Requirements & reproducibility

- **R ≥ 4.0** (tested on 4.5 / 4.6), base R only; `biomaRt`/`jsonlite` optional.
- **`curl`** required for all three downloads.
- Outbound **HTTPS** required; ~2 GB free disk recommended.
- Download date, host/URL, and row counts are recorded in
  `data/raw/download_metadata.rds` and surfaced in the generated report.
- For a frozen snapshot, archive `data/raw/` after a successful download.

See `requirements_notes.md` for troubleshooting (BioMart attribute limits, mirror timeouts,
HGNC 404 fallback, R version / biomaRt mismatch).
