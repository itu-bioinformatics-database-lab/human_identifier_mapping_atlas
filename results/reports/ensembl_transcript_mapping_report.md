# Ensembl transcript identifier mappings

Generated: `2026-07-29T19:33:33.037209+00:00`
Source: `data/raw/ensembl_biomart_raw.csv`

## Scope

Adds transcript-level pairwise tables alongside the existing gene / symbol / UniProt atlas:

- `ensembl_transcript_id` → `ensembl_gene_id`
- `ensembl_transcript_id` → `uniprot_id`
- `ensembl_transcript_id` → `gene_symbol`

Schema matches other pairwise JSON files (`source_id`, `target_id`, `source_id_type`, `target_id_type`, `mapping_source`, `source_database`).

## Atlas coverage

- Unique Ensembl transcripts observed: **670,670**
- With ≥1 Ensembl gene: **670,670**
- With ≥1 gene symbol: **558,398**
- With ≥1 UniProt accession: **133,934**

## Pairwise stats

| source | target | n_rows | n_source | n_target |
|--------|--------|-------:|---------:|---------:|
| ensembl_transcript_id | ensembl_gene_id | 670,670 | 670,670 | 86,411 |
| ensembl_transcript_id | uniprot_id | 141,653 | 133,934 | 86,108 |
| ensembl_transcript_id | gene_symbol | 558,533 | 558,398 | 41,860 |

## Mayo transcript smoke coverage

File: `/home/bioinformatics/Desktop/ad-multiomics-graph/benchmark-data/counts/Mayo_TCX_gene_tx_protein/transcript.csv`

- ENST features: **208,244**
- → Ensembl gene: **198,637** (95.39%)
- → gene symbol: **180,638** (86.74%)
- → UniProt: **94,171** (45.22%)
- Unmapped to gene: **9,607**

