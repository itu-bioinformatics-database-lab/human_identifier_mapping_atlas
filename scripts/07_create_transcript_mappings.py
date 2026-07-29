#!/usr/bin/env python3
"""Build Ensembl transcript pairwise mappings from BioMart raw.

Reads ``data/raw/ensembl_biomart_raw.csv`` (already downloaded by step 01) and
writes list-of-record JSON edges in the same schema as step 03:

  - ensembl_transcript_to_ensembl_gene.json
  - ensembl_transcript_to_uniprot.json
  - ensembl_transcript_to_gene_symbol.json

Also writes transcript mapping stats / unmapped summaries under
``results/stats/`` and a tracked report under ``results/reports/``.

Does not modify gene / symbol / UniProt tables produced by the R pipeline.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "ensembl_biomart_raw.csv"
OUT_TABLES = ROOT / "results" / "mapping_tables"
OUT_STATS = ROOT / "results" / "stats"
OUT_UNMAPPED = ROOT / "data" / "unmapped"
OUT_REPORTS = ROOT / "results" / "reports"

MAPPING_SOURCE = "ensembl_biomart"
SOURCE_DATABASE = "ensembl_biomart"

_ENST_RE = re.compile(r"^ENST\d{11}(?:\.\d+)?$", re.IGNORECASE)
_ENSG_RE = re.compile(r"^ENSG\d{11}(?:\.\d+)?$", re.IGNORECASE)


def blank(x) -> str | None:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "na", "-", "null"}:
        return None
    return s


def strip_version(x: str) -> str:
    if "." in x and x.upper().startswith("ENS"):
        return x.split(".", 1)[0]
    return x


def split_multi(x: str | None) -> list[str]:
    v = blank(x)
    if v is None:
        return []
    parts = re.split(r"[,|;/\s]+", v)
    out = []
    for p in parts:
        p = blank(p)
        if p:
            out.append(p)
    return out


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n")


def write_edges(path: Path, edges: list[dict]) -> None:
    # Compact JSON (no pretty indent) — tables are large.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(edges, f, separators=(",", ":"))
        f.write("\n")


def pairwise_stats(edges: list[dict], source_type: str, target_type: str) -> dict:
    src_to_tgt: dict[str, set[str]] = defaultdict(set)
    tgt_to_src: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        s, t = e["source_id"], e["target_id"]
        src_to_tgt[s].add(t)
        tgt_to_src[t].add(s)
    n_one_to_one = n_one_to_many = n_many_to_one = n_many_to_many = 0
    for s, tgts in src_to_tgt.items():
        many_t = len(tgts) > 1
        many_s = any(len(tgt_to_src[t]) > 1 for t in tgts)
        if not many_t and not many_s:
            n_one_to_one += 1
        elif many_t and not many_s:
            n_one_to_many += 1
        elif not many_t and many_s:
            n_many_to_one += 1
        else:
            n_many_to_many += 1
    return {
        "source_id_type": source_type,
        "target_id_type": target_type,
        "n_unique_source_ids": len(src_to_tgt),
        "n_unique_target_ids": len(tgt_to_src),
        "n_mapping_rows": len(edges),
        "n_one_to_one": n_one_to_one,
        "n_one_to_many": n_one_to_many,
        "n_many_to_one": n_many_to_one,
        "n_many_to_many": n_many_to_many,
        "mapping_source": MAPPING_SOURCE,
        "source_database": SOURCE_DATABASE,
    }


def build_maps(raw_path: Path):
    tx2gene: dict[str, set[str]] = defaultdict(set)
    tx2sym: dict[str, set[str]] = defaultdict(set)
    tx2up: dict[str, set[str]] = defaultdict(set)
    all_tx: set[str] = set()

    with raw_path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "external_gene_name",
            "hgnc_symbol",
            "uniprotswissprot",
            "uniprotsptrembl",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"BioMart raw missing columns: {sorted(missing)}")

        for row in reader:
            tx = blank(row.get("ensembl_transcript_id"))
            if tx is None:
                continue
            tx = strip_version(tx)
            if not _ENST_RE.match(tx):
                continue
            all_tx.add(tx)

            gene = blank(row.get("ensembl_gene_id"))
            if gene is not None:
                gene = strip_version(gene)
                if _ENSG_RE.match(gene):
                    tx2gene[tx].add(gene)

            sym = blank(row.get("hgnc_symbol")) or blank(row.get("external_gene_name"))
            if sym is not None:
                tx2sym[tx].add(sym)

            # Prefer Swiss-Prot; still keep TrEMBL so one-to-many is preserved.
            for col in ("uniprotswissprot", "uniprotsptrembl"):
                for uid in split_multi(row.get(col)):
                    tx2up[tx].add(uid)

    return all_tx, tx2gene, tx2sym, tx2up


def edges_from_map(
    mapping: dict[str, set[str]],
    source_type: str,
    target_type: str,
) -> list[dict]:
    edges: list[dict] = []
    for sid in sorted(mapping):
        for tid in sorted(mapping[sid]):
            edges.append(
                {
                    "source_id": sid,
                    "target_id": tid,
                    "source_id_type": source_type,
                    "target_id_type": target_type,
                    "mapping_source": MAPPING_SOURCE,
                    "source_database": SOURCE_DATABASE,
                }
            )
    return edges


def mayo_coverage(
    mayo_header_csv: Path | None,
    tx2gene: dict[str, set[str]],
    tx2sym: dict[str, set[str]],
    tx2up: dict[str, set[str]],
) -> dict | None:
    if mayo_header_csv is None or not mayo_header_csv.is_file():
        return None
    with mayo_header_csv.open(newline="") as f:
        cols = next(csv.reader(f))
    feats = [c for c in cols if c not in {"sample_id", "diagnosis", "Sample ID"}]
    ids = []
    for c in feats:
        c = strip_version(c.strip())
        if _ENST_RE.match(c):
            ids.append(c)
    n = len(ids)
    if n == 0:
        return {"n_features": 0}
    g = sum(1 for i in ids if i in tx2gene)
    s = sum(1 for i in ids if i in tx2sym)
    u = sum(1 for i in ids if i in tx2up)
    unmapped = [i for i in ids if i not in tx2gene]
    return {
        "source_file": str(mayo_header_csv),
        "n_ensembl_transcript_features": n,
        "mapped_to_ensembl_gene": g,
        "mapped_to_gene_symbol": s,
        "mapped_to_uniprot": u,
        "unmapped_to_ensembl_gene": n - g,
        "mapping_rate_to_ensembl_gene": round(g / n, 6),
        "mapping_rate_to_gene_symbol": round(s / n, 6),
        "mapping_rate_to_uniprot": round(u / n, 6),
        "unmapped_to_ensembl_gene_examples": unmapped[:25],
    }


def write_report(
    stats_rows: list[dict],
    all_tx: set[str],
    tx2gene: dict[str, set[str]],
    tx2sym: dict[str, set[str]],
    tx2up: dict[str, set[str]],
    mayo: dict | None,
) -> Path:
    lines = [
        "# Ensembl transcript identifier mappings",
        "",
        f"Generated: `{datetime.now(timezone.utc).isoformat()}`",
        f"Source: `{RAW.relative_to(ROOT)}`",
        "",
        "## Scope",
        "",
        "Adds transcript-level pairwise tables alongside the existing gene / symbol / UniProt atlas:",
        "",
        "- `ensembl_transcript_id` → `ensembl_gene_id`",
        "- `ensembl_transcript_id` → `uniprot_id`",
        "- `ensembl_transcript_id` → `gene_symbol`",
        "",
        "Schema matches other pairwise JSON files (`source_id`, `target_id`, "
        "`source_id_type`, `target_id_type`, `mapping_source`, `source_database`).",
        "",
        "## Atlas coverage",
        "",
        f"- Unique Ensembl transcripts observed: **{len(all_tx):,}**",
        f"- With ≥1 Ensembl gene: **{len(tx2gene):,}**",
        f"- With ≥1 gene symbol: **{len(tx2sym):,}**",
        f"- With ≥1 UniProt accession: **{len(tx2up):,}**",
        "",
        "## Pairwise stats",
        "",
        "| source | target | n_rows | n_source | n_target |",
        "|--------|--------|-------:|---------:|---------:|",
    ]
    for st in stats_rows:
        lines.append(
            f"| {st['source_id_type']} | {st['target_id_type']} | "
            f"{st['n_mapping_rows']:,} | {st['n_unique_source_ids']:,} | "
            f"{st['n_unique_target_ids']:,} |"
        )
    lines += ["", "## Mayo transcript smoke coverage", ""]
    if mayo and mayo.get("n_ensembl_transcript_features"):
        lines += [
            f"File: `{mayo['source_file']}`",
            "",
            f"- ENST features: **{mayo['n_ensembl_transcript_features']:,}**",
            f"- → Ensembl gene: **{mayo['mapped_to_ensembl_gene']:,}** "
            f"({100 * mayo['mapping_rate_to_ensembl_gene']:.2f}%)",
            f"- → gene symbol: **{mayo['mapped_to_gene_symbol']:,}** "
            f"({100 * mayo['mapping_rate_to_gene_symbol']:.2f}%)",
            f"- → UniProt: **{mayo['mapped_to_uniprot']:,}** "
            f"({100 * mayo['mapping_rate_to_uniprot']:.2f}%)",
            f"- Unmapped to gene: **{mayo['unmapped_to_ensembl_gene']:,}**",
            "",
        ]
    else:
        lines += ["Mayo transcript header CSV not provided; skipped.", ""]
    path = OUT_REPORTS / "ensembl_transcript_mapping_report.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mayo-transcript-csv",
        type=Path,
        default=Path(
            "/home/bioinformatics/Desktop/ad-multiomics-graph/benchmark-data/"
            "counts/Mayo_TCX_gene_tx_protein/transcript.csv"
        ),
        help="Optional Mayo transcript matrix (header used for coverage).",
    )
    p.add_argument(
        "--raw",
        type=Path,
        default=RAW,
        help="Path to ensembl_biomart_raw.csv",
    )
    args = p.parse_args()

    if not args.raw.is_file():
        raise SystemExit(f"Missing BioMart raw file: {args.raw}")

    print(f"Reading {args.raw} ...", flush=True)
    all_tx, tx2gene, tx2sym, tx2up = build_maps(args.raw)
    print(
        f"  transcripts={len(all_tx):,}  "
        f"with_gene={len(tx2gene):,}  with_symbol={len(tx2sym):,}  with_uniprot={len(tx2up):,}",
        flush=True,
    )

    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_STATS.mkdir(parents=True, exist_ok=True)
    OUT_UNMAPPED.mkdir(parents=True, exist_ok=True)

    specs = [
        (tx2gene, "ensembl_transcript_id", "ensembl_gene_id", "ensembl_transcript_to_ensembl_gene.json"),
        (tx2up, "ensembl_transcript_id", "uniprot_id", "ensembl_transcript_to_uniprot.json"),
        (tx2sym, "ensembl_transcript_id", "gene_symbol", "ensembl_transcript_to_gene_symbol.json"),
    ]
    stats_rows = []
    for mapping, stype, ttype, fname in specs:
        edges = edges_from_map(mapping, stype, ttype)
        out = OUT_TABLES / fname
        print(f"  writing {out.name}: {len(edges):,} edges ...", flush=True)
        write_edges(out, edges)
        stats_rows.append(pairwise_stats(edges, stype, ttype))

    write_json(OUT_STATS / "ensembl_transcript_pairwise_stats.json", stats_rows)

    # Unmapped transcripts relative to each target type
    unmapped_specs = [
        (all_tx - set(tx2gene), "ensembl_transcripts_not_mapped_to_ensembl_gene.txt"),
        (all_tx - set(tx2up), "ensembl_transcripts_not_mapped_to_uniprot.txt"),
        (all_tx - set(tx2sym), "ensembl_transcripts_not_mapped_to_gene_symbol.txt"),
    ]
    unmapped_summary = []
    for missing, fname in unmapped_specs:
        path = OUT_UNMAPPED / fname
        path.write_text("\n".join(sorted(missing)) + ("\n" if missing else ""))
        unmapped_summary.append(
            {
                "file": fname,
                "n_unmapped": len(missing),
                "n_total_transcripts": len(all_tx),
            }
        )
    write_json(OUT_STATS / "ensembl_transcript_unmapped_summary.json", unmapped_summary)

    mayo = mayo_coverage(args.mayo_transcript_csv, tx2gene, tx2sym, tx2up)
    if mayo:
        write_json(OUT_STATS / "mayo_transcript_mapping_coverage.json", mayo)
        print(
            f"  Mayo coverage: n={mayo.get('n_ensembl_transcript_features')} "
            f"gene={mayo.get('mapped_to_ensembl_gene')} "
            f"symbol={mayo.get('mapped_to_gene_symbol')} "
            f"uniprot={mayo.get('mapped_to_uniprot')}",
            flush=True,
        )

    report = write_report(stats_rows, all_tx, tx2gene, tx2sym, tx2up, mayo)
    print(f"Wrote report {report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
