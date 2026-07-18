# human_identifier_mapping_atlas

Focused human gene/protein identifier mapping atlas for three ID types only:

- Ensembl gene IDs
- Gene symbols
- UniProt accessions

All valid one-to-many, many-to-one, and many-to-many mappings are preserved.
Other ID types are ignored. Mapping tables are written as **JSON**.

## How to run

Assumes raw downloads already exist under `data/raw/`.

```bash
cd ~/Desktop/human_identifier_mapping_atlas
Rscript scripts/run_all.R
```

Or step by step:

```bash
Rscript scripts/02_process_sources.R
Rscript scripts/03_create_pairwise_mappings.R
Rscript scripts/04_generate_unmapped_lists.R
Rscript scripts/05_generate_report.R
```

To (re)download sources:

```bash
FORCE_REDOWNLOAD=1 Rscript scripts/01_download_sources.R
```

## Outputs

| Path | Description |
|------|-------------|
| `results/mapping_tables/*.json` | Pairwise + master mappings |
| `results/stats/` | Pairwise stats, ID counts, unmapped summary |
| `data/unmapped/*.txt` | One ID per line |
| `results/reports/identifier_mapping_report.md` | Summary report |

### Pairwise JSON files

- `ensembl_gene_to_gene_symbol.json`
- `ensembl_gene_to_uniprot.json`
- `gene_symbol_to_ensembl_gene.json`
- `gene_symbol_to_uniprot.json`
- `uniprot_to_ensembl_gene.json`
- `uniprot_to_gene_symbol.json`
- `master_identifier_mappings_long.json`
- `master_identifier_annotations_wide.json`

## Artifact mappers (UniProt / HGNC downloads)

If you place source downloads under `artifacts/mappers/`:

```text
artifacts/mappers/UniProt/HUMAN_9606_idmapping.dat.gz
artifacts/mappers/HGNC/hgnc_complete_set.txt
```

Build BioMart-style focused mappings (Ensembl gene / gene symbol / UniProt) in
**separate folders** for each source:

```bash
python3 scripts/06_process_artifact_mappers.py
```

This writes, under each source folder:

- `mappings/*.json`
- `unmapped/*.txt`
- `stats/*`
- `reports/identifier_mapping_report.md`

Raw download files are never deleted.

## Data sources & large files

How each database is downloaded (Ensembl BioMart, HGNC, UniProt), where files land,
and what is excluded from git is documented in **[DATA_SOURCES.md](DATA_SOURCES.md)**.

Raw downloads and generated mapping tables are **gitignored** (often >100 MB). After cloning:

```bash
Rscript scripts/01_download_sources.R
Rscript scripts/run_all.R
```

## Notes

- If UniProt proteome download failed, UniProt IDs are still taken from Ensembl BioMart and HGNC cross-references.
- Raw downloads under `data/raw/` are never cleared by processing scripts.
