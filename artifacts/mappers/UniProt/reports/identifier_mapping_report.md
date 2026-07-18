# UniProt Identifier Mappings — Report

**Report generated:** 2026-07-08 16:05:36

## Scope

Focused BioMart-style mappings among:

- `ensembl_gene_id`
- `gene_symbol`
- `uniprot_id`

One-to-many / many-to-one / many-to-many mappings are preserved.

## Unique ID counts

| ID type | Count |
|---------|------:|
| ensembl_gene_id | 23,740 |
| gene_symbol | 25,496 |
| uniprot_id | 210,709 |

## Pairwise mapping counts

| Source | Target | Unique sources | Unique targets | Mapping rows | 1:1 | 1:many | many:1 | many:many |
|--------|--------|---------------:|---------------:|-------------:|----:|-------:|-------:|----------:|
| ensembl_gene_id | gene_symbol | 23,164 | 19,944 | 24,333 | 17,712 | 146 | 4,921 | 1,554 |
| ensembl_gene_id | uniprot_id | 23,740 | 146,520 | 153,608 | 4,249 | 137,704 | 1,998 | 9,657 |
| gene_symbol | ensembl_gene_id | 19,944 | 23,164 | 24,333 | 17,712 | 4,921 | 146 | 1,554 |
| gene_symbol | uniprot_id | 25,496 | 171,710 | 171,835 | 8,442 | 163,200 | 135 | 58 |
| uniprot_id | ensembl_gene_id | 146,520 | 23,740 | 153,608 | 4,249 | 1,998 | 137,704 | 9,657 |
| uniprot_id | gene_symbol | 171,710 | 25,496 | 171,835 | 8,442 | 135 | 163,200 | 58 |

## Mapping rates / unmapped

| Source | Target | Total | Mapped | Unmapped | Mapping rate |
|--------|--------|------:|-------:|---------:|-------------:|
| ensembl_gene_id | gene_symbol | 23,740 | 23,164 | 576 | 97.57% |
| ensembl_gene_id | uniprot_id | 23,740 | 23,740 | 0 | 100.00% |
| gene_symbol | ensembl_gene_id | 25,496 | 19,944 | 5,552 | 78.22% |
| gene_symbol | uniprot_id | 25,496 | 25,496 | 0 | 100.00% |
| uniprot_id | ensembl_gene_id | 210,709 | 146,520 | 64,189 | 69.54% |
| uniprot_id | gene_symbol | 210,709 | 171,710 | 38,999 | 81.49% |

## Notes

- Source file: HUMAN_9606_idmapping.dat.gz
- Gene symbols taken from UniProt `Gene_Name` (primary), not Gene_Synonym.
- Ensembl IDs taken from UniProt `Ensembl` and version suffixes stripped.
- Isoform accessions (e.g. P31946-1) collapsed to canonical accession (P31946).

## Output locations

- mappings: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/UniProt/mappings`
- unmapped: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/UniProt/unmapped`
- stats: `/home/bioinformatics/Desktop/human_identifier_mapping_atlas/artifacts/mappers/UniProt/stats`
