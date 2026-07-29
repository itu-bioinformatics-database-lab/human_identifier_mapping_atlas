#!/usr/bin/env python3
"""Validate Ensembl transcript mapping tables produced by step 07."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "mapping_tables"
STATS = ROOT / "results" / "stats"

REQUIRED = {
    "ensembl_transcript_to_ensembl_gene.json": ("ensembl_transcript_id", "ensembl_gene_id"),
    "ensembl_transcript_to_uniprot.json": ("ensembl_transcript_id", "uniprot_id"),
    "ensembl_transcript_to_gene_symbol.json": ("ensembl_transcript_id", "gene_symbol"),
}


class TranscriptMappingTests(unittest.TestCase):
    def test_tables_exist_and_schema(self):
        for fname, (stype, ttype) in REQUIRED.items():
            path = TABLES / fname
            self.assertTrue(path.is_file(), f"missing {path}")
            records = json.loads(path.read_text())
            self.assertIsInstance(records, list)
            self.assertGreater(len(records), 0, fname)
            row = records[0]
            for key in (
                "source_id",
                "target_id",
                "source_id_type",
                "target_id_type",
                "mapping_source",
                "source_database",
            ):
                self.assertIn(key, row)
            self.assertEqual(row["source_id_type"], stype)
            self.assertEqual(row["target_id_type"], ttype)
            self.assertTrue(str(row["source_id"]).startswith("ENST"))

    def test_known_mitochondrial_transcript(self):
        # Compact JSON from step 07 — spot-check a known BioMart row.
        raw = (TABLES / "ensembl_transcript_to_ensembl_gene.json").read_text()
        self.assertIn('"source_id":"ENST00000361390","target_id":"ENSG00000198888"', raw)

    def test_mayo_coverage_stats(self):
        path = STATS / "mayo_transcript_mapping_coverage.json"
        self.assertTrue(path.is_file(), "run 07 with --mayo-transcript-csv first")
        cov = json.loads(path.read_text())
        n = cov["n_ensembl_transcript_features"]
        self.assertGreater(n, 100_000)
        self.assertGreaterEqual(cov["mapping_rate_to_ensembl_gene"], 0.90)
        self.assertGreaterEqual(cov["mapped_to_ensembl_gene"], int(0.90 * n))
        # UniProt is expected to be partial (many non-coding / no Swiss-Prot).
        self.assertGreater(cov["mapped_to_uniprot"], 0)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TranscriptMappingTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
