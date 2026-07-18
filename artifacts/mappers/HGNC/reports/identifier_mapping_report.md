# HGNC Identifier Mappings — Report

**Report generated:** 2026-07-08 16:05:38

## Scope

Focused BioMart-style mappings among:

- `ensembl_gene_id`
- `gene_symbol`
- `uniprot_id`

One-to-many / many-to-one / many-to-many mappings are preserved.

## Unique ID counts

| ID type | Count |
|---------|------:|
| ensembl_gene_id | 42,334 |
| gene_symbol | 45,019 |
| uniprot_id | 20,251 |

## Pairwise mapping counts

| Source | Target | Unique sources | Unique targets | Mapping rows | 1:1 | 1:many | many:1 | many:many |
|--------|--------|---------------:|---------------:|-------------:|----:|-------:|-------:|----------:|
| ensembl_gene_id | gene_symbol | 42,334 | 42,337 | 42,337 | 42,331 | 6 | 0 | 0 |
| ensembl_gene_id | uniprot_id | 20,200 | 20,127 | 20,251 | 19,965 | 93 | 191 | 2 |
| gene_symbol | ensembl_gene_id | 42,337 | 42,334 | 42,337 | 42,331 | 0 | 6 | 0 |
| gene_symbol | uniprot_id | 20,288 | 20,251 | 20,377 | 20,036 | 144 | 195 | 2 |
| uniprot_id | ensembl_gene_id | 20,127 | 20,200 | 20,251 | 19,965 | 191 | 93 | 2 |
| uniprot_id | gene_symbol | 20,251 | 20,288 | 20,377 | 20,036 | 195 | 144 | 2 |

## Mapping rates / unmapped

| Source | Target | Total | Mapped | Unmapped | Mapping rate |
|--------|--------|------:|-------:|---------:|-------------:|
| ensembl_gene_id | gene_symbol | 42,334 | 42,334 | 0 | 100.00% |
| ensembl_gene_id | uniprot_id | 42,334 | 20,200 | 22,134 | 47.72% |
| gene_symbol | ensembl_gene_id | 45,019 | 42,337 | 2,682 | 94.04% |
| gene_symbol | uniprot_id | 45,019 | 20,288 | 24,731 | 45.07% |
| uniprot_id | ensembl_gene_id | 20,251 | 20,127 | 124 | 99.39% |
| uniprot_id | gene_symbol | 20,251 | 20,251 | 0 | 100.00% |

## Notes

- Source file: hgnc_complete_set.txt
- Gene symbols are HGNC approved symbols.
- Multi-value UniProt fields are split into long-format mappings.
- Aliases / previous symbols are ignored for this focused three-ID atlas.

## Output locations

- mappings: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/HGNC/mappings`
- unmapped: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/HGNC/unmapped`
- stats: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/HGNC/stats`
