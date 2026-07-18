# Requirements notes — human_identifier_mapping_atlas

## R version

Tested with R ≥ 4.0. Works on R 4.5 and R 4.6.

## Dependencies

### Required system tools

- `curl` (BioMart XML POST + UniProt stream download)

### R packages

**No CRAN packages are required.** The pipeline is implemented in base R with a
project-root path helper (`here_path()` in `scripts/00_utils.R`).

Optional:

| Package | Purpose |
|---------|---------|
| `biomaRt` (Bioconductor) | Faster Ensembl downloads when it loads successfully |
| `here`, `httr` | Not required; historical / optional |

If `biomaRt` is unavailable or fails to load (e.g. R 4.6 vs packages compiled
for R 4.5), step 01 automatically falls back to the BioMart HTTP XML API.

## Data sources

### 1. Ensembl BioMart

- Dataset: `hsapiens_gene_ensembl`
- Hosts tried: `www.ensembl.org`, `useast.ensembl.org`, `asia.ensembl.org`
- Attributes are fetched in **three queries** (core / UniProt+Entrez / RefSeq)
  and merged, to avoid Ensembl’s “Too many attributes selected for External
  References” error.

### 2. HGNC

- Primary URL: `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
- Fallback: genenames.org custom download endpoint
- Format: TSV

### 3. UniProt

- Endpoint: `https://rest.uniprot.org/uniprotkb/stream` (TSV)
- Query: `organism_id:9606`
- Includes reviewed (Swiss-Prot) and unreviewed (TrEMBL)

## Network requirements

All three sources require outbound HTTPS. If downloads fail:

1. Check firewall/proxy settings
2. Re-run `FORCE_REDOWNLOAD=1 Rscript scripts/01_download_sources.R`
3. Inspect `data/raw/download_metadata.rds`

## Disk space

| File | Approximate size |
|------|-----------------|
| `ensembl_biomart_raw.csv` | 50–150 MB |
| `hgnc_complete_raw.csv` | 5–10 MB |
| `uniprot_human_raw.csv` | 100–300 MB |
| Master mapping tables | 100–500 MB |

Ensure ≥ 2 GB free disk space.

## Runtime estimates

| Step | Typical duration |
|------|-----------------|
| 01 Download | 10–40 min |
| 02 Process | 2–10 min |
| 03 Pairwise mappings | 3–15 min |
| 04 Unmapped lists | 1–5 min |
| 05 Report | < 1 min |

## Reproducibility

- Download date is recorded in `data/raw/download_metadata.rds` and the markdown report
- For frozen snapshots, archive `data/raw/` after a successful download

## Troubleshooting

| Issue | Solution |
|-------|----------|
| BioMart "too many attributes" | Handled by split queries in `01_download_sources.R` |
| BioMart timeout | Re-run; script tries alternate Ensembl mirrors |
| UniProt HTTP/timeout | Stream endpoint is default; check connection stability |
| HGNC 404 | Script falls back to an alternate genenames.org URL |
| R 4.6 cannot load R 4.5 `biomaRt` | HTTP fallback is used automatically |
| Out of memory | Use a machine with ≥ 8 GB RAM |
