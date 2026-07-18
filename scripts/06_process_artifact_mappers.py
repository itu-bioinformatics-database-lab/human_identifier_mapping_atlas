#!/usr/bin/env python3
"""Build BioMart-style Ensembl / gene_symbol / UniProt mappings from
artifacts/mappers/UniProt and artifacts/mappers/HGNC.

Outputs (per source, under that source folder):
  mappings/*.json
  unmapped/*.txt
  stats/*.csv + *.json
  reports/identifier_mapping_report.md

Raw download files are never deleted.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts" / "mappers"

FOCUS_PAIRS = [
    ("ensembl_gene_id", "gene_symbol"),
    ("ensembl_gene_id", "uniprot_id"),
    ("gene_symbol", "ensembl_gene_id"),
    ("gene_symbol", "uniprot_id"),
    ("uniprot_id", "ensembl_gene_id"),
    ("uniprot_id", "gene_symbol"),
]

PAIR_FILENAMES = {
    ("ensembl_gene_id", "gene_symbol"): "ensembl_gene_to_gene_symbol.json",
    ("ensembl_gene_id", "uniprot_id"): "ensembl_gene_to_uniprot.json",
    ("gene_symbol", "ensembl_gene_id"): "gene_symbol_to_ensembl_gene.json",
    ("gene_symbol", "uniprot_id"): "gene_symbol_to_uniprot.json",
    ("uniprot_id", "ensembl_gene_id"): "uniprot_to_ensembl_gene.json",
    ("uniprot_id", "gene_symbol"): "uniprot_to_gene_symbol.json",
}

UNMAPPED_SPECS = [
    ("ensembl_gene_id", "gene_symbol", "ensembl_genes_not_mapped_to_gene_symbol.txt"),
    ("ensembl_gene_id", "uniprot_id", "ensembl_genes_not_mapped_to_uniprot.txt"),
    ("gene_symbol", "ensembl_gene_id", "gene_symbols_not_mapped_to_ensembl_gene.txt"),
    ("gene_symbol", "uniprot_id", "gene_symbols_not_mapped_to_uniprot.txt"),
    ("uniprot_id", "ensembl_gene_id", "uniprots_not_mapped_to_ensembl_gene.txt"),
    ("uniprot_id", "gene_symbol", "uniprots_not_mapped_to_gene_symbol.txt"),
]


def ensure_dirs(base: Path) -> dict[str, Path]:
    dirs = {
        "mappings": base / "mappings",
        "unmapped": base / "unmapped",
        "stats": base / "stats",
        "reports": base / "reports",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def blank_to_none(x: str | None) -> str | None:
    if x is None:
        return None
    x = str(x).strip()
    if not x or x in {"-", "NA", "None", "null"}:
        return None
    return x


def split_multi(x: str | None) -> list[str]:
    v = blank_to_none(x)
    if v is None:
        return []
    parts = re.split(r"[,|;/\s]+", v)
    out = []
    seen = set()
    for p in parts:
        p = p.strip()
        if not p or p in {"-", "NA"}:
            continue
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def strip_ensembl_version(x: str) -> str:
    # ENSG00000166913.14 -> ENSG00000166913
    return re.sub(r"\.\d+$", "", x)


def canonical_uniprot(acc: str) -> str:
    # P31946-1 -> P31946
    return acc.split("-", 1)[0]


def write_json(path: Path, obj, pretty: bool = False) -> None:
    with path.open("w", encoding="utf-8") as f:
        if pretty:
            json.dump(obj, f, indent=2, ensure_ascii=False)
        else:
            json.dump(obj, f, separators=(",", ":"), ensure_ascii=False)
    n = len(obj) if isinstance(obj, list) else None
    print(f"  Wrote {path} ({n} records)" if n is not None else f"  Wrote {path}")


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  Wrote {path} ({len(rows)} rows)")


def edges_from_pairs(
    pairs: set[tuple[str, str]],
    source_type: str,
    target_type: str,
    mapping_source: str,
    source_database: str,
) -> list[dict]:
    return [
        {
            "source_id": s,
            "target_id": t,
            "source_id_type": source_type,
            "target_id_type": target_type,
            "mapping_source": mapping_source,
            "source_database": source_database,
        }
        for s, t in sorted(pairs)
    ]


def compute_pairwise_stats(edges: list[dict], source_type: str, target_type: str) -> dict:
    if not edges:
        return {
            "source_id_type": source_type,
            "target_id_type": target_type,
            "n_unique_source_ids": 0,
            "n_unique_target_ids": 0,
            "n_mapping_rows": 0,
            "n_one_to_one": 0,
            "n_one_to_many": 0,
            "n_many_to_one": 0,
            "n_many_to_many": 0,
            "mapping_rate": None,
            "source_database_counts": "",
        }

    src_to_tgt: dict[str, set[str]] = defaultdict(set)
    tgt_to_src: dict[str, set[str]] = defaultdict(set)
    db_counts: dict[str, int] = defaultdict(int)
    for e in edges:
        src_to_tgt[e["source_id"]].add(e["target_id"])
        tgt_to_src[e["target_id"]].add(e["source_id"])
        db_counts[e.get("source_database", "")] += 1

    n_one_to_one = n_one_to_many = n_many_to_one = n_many_to_many = 0
    for e in edges:
        sc = len(src_to_tgt[e["source_id"]])
        tc = len(tgt_to_src[e["target_id"]])
        if sc == 1 and tc == 1:
            n_one_to_one += 1
        elif sc > 1 and tc == 1:
            n_one_to_many += 1
        elif sc == 1 and tc > 1:
            n_many_to_one += 1
        else:
            n_many_to_many += 1

    n_src = len(src_to_tgt)
    return {
        "source_id_type": source_type,
        "target_id_type": target_type,
        "n_unique_source_ids": n_src,
        "n_unique_target_ids": len(tgt_to_src),
        "n_mapping_rows": len(edges),
        "n_one_to_one": n_one_to_one,
        "n_one_to_many": n_one_to_many,
        "n_many_to_one": n_many_to_one,
        "n_many_to_many": n_many_to_many,
        "mapping_rate": round(len(edges) / max(n_src, 1), 4),
        "source_database_counts": ";".join(f"{k}:{v}" for k, v in sorted(db_counts.items())),
    }


def write_outputs(
    out_dirs: dict[str, Path],
    edge_store: dict[tuple[str, str], list[dict]],
    id_pools: dict[str, set[str]],
    source_label: str,
    extra_notes: list[str] | None = None,
) -> None:
    # Pairwise JSON
    all_edges: list[dict] = []
    stats_rows: list[dict] = []
    for src, tgt in FOCUS_PAIRS:
        edges = edge_store.get((src, tgt), [])
        fname = PAIR_FILENAMES[(src, tgt)]
        write_json(out_dirs["mappings"] / fname, edges)
        all_edges.extend(edges)
        stats_rows.append(compute_pairwise_stats(edges, src, tgt))

    # Master long / wide
    write_json(out_dirs["mappings"] / "master_identifier_mappings_long.json", all_edges)

    # Wide annotations: rebuild unique triples from available pools/edges
    wide_rows = set()
    for (src_t, tgt_t), edges in edge_store.items():
        for e in edges:
            row = {"ensembl_gene_id": None, "gene_symbol": None, "uniprot_id": None, "source_database": e["source_database"]}
            row[src_t] = e["source_id"]
            row[tgt_t] = e["target_id"]
            # Only keep rows that have at least two of the three IDs filled via this edge
            wide_rows.add((row["ensembl_gene_id"], row["gene_symbol"], row["uniprot_id"], row["source_database"]))

    wide_list = [
        {
            "ensembl_gene_id": a,
            "gene_symbol": b,
            "uniprot_id": c,
            "source_database": d,
        }
        for a, b, c, d in sorted(wide_rows, key=lambda x: (x[0] or "", x[1] or "", x[2] or "", x[3] or ""))
    ]
    write_json(out_dirs["mappings"] / "master_identifier_annotations_wide.json", wide_list)

    # Stats
    id_counts = [
        {"id_type": k, "count": len(v)}
        for k, v in [
            ("ensembl_gene_id", id_pools["ensembl_gene_id"]),
            ("gene_symbol", id_pools["gene_symbol"]),
            ("uniprot_id", id_pools["uniprot_id"]),
        ]
    ]
    write_csv(out_dirs["stats"] / "id_type_counts.csv", id_counts, ["id_type", "count"])
    write_json(out_dirs["stats"] / "id_type_counts.json", id_counts, pretty=True)
    write_csv(
        out_dirs["stats"] / "pairwise_mapping_stats.csv",
        stats_rows,
        list(stats_rows[0].keys()),
    )
    write_json(out_dirs["stats"] / "pairwise_mapping_stats.json", stats_rows, pretty=True)

    # Unmapped
    mapped_index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for e in all_edges:
        mapped_index[(e["source_id_type"], e["target_id_type"])].add(e["source_id"])

    unmapped_rows = []
    for src, tgt, fname in UNMAPPED_SPECS:
        pool = id_pools[src]
        mapped = mapped_index.get((src, tgt), set())
        unmapped = sorted(pool - mapped)
        path = out_dirs["unmapped"] / fname
        path.write_text("\n".join(unmapped) + ("\n" if unmapped else ""), encoding="utf-8")
        print(f"  Wrote {path} ({len(unmapped)} IDs)")
        total = len(pool)
        n_unmapped = len(unmapped)
        n_mapped = total - n_unmapped
        unmapped_rows.append(
            {
                "source_id_type": src,
                "target_id_type": tgt,
                "total_source_ids": total,
                "mapped_source_ids": n_mapped,
                "unmapped_source_ids": n_unmapped,
                "mapping_rate": round(n_mapped / total, 4) if total else None,
                "unmapped_rate": round(n_unmapped / total, 4) if total else None,
                "unmapped_file": fname,
            }
        )
    write_csv(out_dirs["stats"] / "unmapped_summary.csv", unmapped_rows, list(unmapped_rows[0].keys()))
    write_json(out_dirs["stats"] / "unmapped_summary.json", unmapped_rows, pretty=True)

    # Report
    lines = [
        f"# {source_label} Identifier Mappings — Report",
        "",
        f"**Report generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S %Z')}".rstrip(),
        "",
        "## Scope",
        "",
        "Focused BioMart-style mappings among:",
        "",
        "- `ensembl_gene_id`",
        "- `gene_symbol`",
        "- `uniprot_id`",
        "",
        "One-to-many / many-to-one / many-to-many mappings are preserved.",
        "",
        "## Unique ID counts",
        "",
        "| ID type | Count |",
        "|---------|------:|",
    ]
    for r in id_counts:
        lines.append(f"| {r['id_type']} | {r['count']:,} |")
    lines += [
        "",
        "## Pairwise mapping counts",
        "",
        "| Source | Target | Unique sources | Unique targets | Mapping rows | 1:1 | 1:many | many:1 | many:many |",
        "|--------|--------|---------------:|---------------:|-------------:|----:|-------:|-------:|----------:|",
    ]
    for r in stats_rows:
        lines.append(
            f"| {r['source_id_type']} | {r['target_id_type']} | "
            f"{r['n_unique_source_ids']:,} | {r['n_unique_target_ids']:,} | {r['n_mapping_rows']:,} | "
            f"{r['n_one_to_one']:,} | {r['n_one_to_many']:,} | {r['n_many_to_one']:,} | {r['n_many_to_many']:,} |"
        )
    lines += [
        "",
        "## Mapping rates / unmapped",
        "",
        "| Source | Target | Total | Mapped | Unmapped | Mapping rate |",
        "|--------|--------|------:|-------:|---------:|-------------:|",
    ]
    for r in unmapped_rows:
        rate = "N/A" if r["mapping_rate"] is None else f"{r['mapping_rate'] * 100:.2f}%"
        lines.append(
            f"| {r['source_id_type']} | {r['target_id_type']} | {r['total_source_ids']:,} | "
            f"{r['mapped_source_ids']:,} | {r['unmapped_source_ids']:,} | {rate} |"
        )
    if extra_notes:
        lines += ["", "## Notes", ""] + [f"- {n}" for n in extra_notes]
    lines += [
        "",
        "## Output locations",
        "",
        f"- mappings: `{out_dirs['mappings']}`",
        f"- unmapped: `{out_dirs['unmapped']}`",
        f"- stats: `{out_dirs['stats']}`",
        "",
    ]
    report = out_dirs["reports"] / "identifier_mapping_report.md"
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Wrote {report}")


def process_uniprot(path: Path) -> None:
    print("\n=== Processing UniProt idmapping ===")
    out_dirs = ensure_dirs(path)

    # Collect per-accession annotations
    gene_names: dict[str, set[str]] = defaultdict(set)
    ensembl_genes: dict[str, set[str]] = defaultdict(set)
    # Also track reverse presence for ID pools
    all_uniprot: set[str] = set()
    all_symbols: set[str] = set()
    all_ensembl: set[str] = set()

    n_lines = 0
    with gzip.open(path / "HUMAN_9606_idmapping.dat.gz", "rt", encoding="utf-8") as f:
        for line in f:
            n_lines += 1
            acc, id_type, value = line.rstrip("\n").split("\t", 2)
            can = canonical_uniprot(acc)
            all_uniprot.add(can)
            if id_type == "Gene_Name":
                v = blank_to_none(value)
                if v:
                    gene_names[can].add(v)
                    all_symbols.add(v)
            elif id_type == "Ensembl":
                # May be ENSG....version
                for part in split_multi(value):
                    eg = strip_ensembl_version(part)
                    if eg.startswith("ENSG"):
                        ensembl_genes[can].add(eg)
                        all_ensembl.add(eg)
            # Gene_Synonym intentionally ignored for primary gene_symbol mapping
            # (keeps approved/primary Gene_Name like BioMart external_gene_name)

    print(f"  Read {n_lines:,} idmapping lines")
    print(f"  Unique UniProt accessions: {len(all_uniprot):,}")
    print(f"  Unique Gene_Name symbols: {len(all_symbols):,}")
    print(f"  Unique Ensembl genes: {len(all_ensembl):,}")

    # Build mapping pair sets
    eg_to_sym: set[tuple[str, str]] = set()
    eg_to_up: set[tuple[str, str]] = set()
    sym_to_eg: set[tuple[str, str]] = set()
    sym_to_up: set[tuple[str, str]] = set()
    up_to_eg: set[tuple[str, str]] = set()
    up_to_sym: set[tuple[str, str]] = set()

    for up in all_uniprot:
        syms = gene_names.get(up, set())
        ensg = ensembl_genes.get(up, set())
        for s in syms:
            up_to_sym.add((up, s))
            sym_to_up.add((s, up))
        for e in ensg:
            up_to_eg.add((up, e))
            eg_to_up.add((e, up))
        for e in ensg:
            for s in syms:
                eg_to_sym.add((e, s))
                sym_to_eg.add((s, e))

    db = "uniprot_idmapping"
    edge_store = {
        ("ensembl_gene_id", "gene_symbol"): edges_from_pairs(eg_to_sym, "ensembl_gene_id", "gene_symbol", db, db),
        ("ensembl_gene_id", "uniprot_id"): edges_from_pairs(eg_to_up, "ensembl_gene_id", "uniprot_id", db, db),
        ("gene_symbol", "ensembl_gene_id"): edges_from_pairs(sym_to_eg, "gene_symbol", "ensembl_gene_id", db, db),
        ("gene_symbol", "uniprot_id"): edges_from_pairs(sym_to_up, "gene_symbol", "uniprot_id", db, db),
        ("uniprot_id", "ensembl_gene_id"): edges_from_pairs(up_to_eg, "uniprot_id", "ensembl_gene_id", db, db),
        ("uniprot_id", "gene_symbol"): edges_from_pairs(up_to_sym, "uniprot_id", "gene_symbol", db, db),
    }
    id_pools = {
        "ensembl_gene_id": all_ensembl,
        "gene_symbol": all_symbols,
        "uniprot_id": all_uniprot,
    }
    write_outputs(
        out_dirs,
        edge_store,
        id_pools,
        source_label="UniProt",
        extra_notes=[
            "Source file: HUMAN_9606_idmapping.dat.gz",
            "Gene symbols taken from UniProt `Gene_Name` (primary), not Gene_Synonym.",
            "Ensembl IDs taken from UniProt `Ensembl` and version suffixes stripped.",
            "Isoform accessions (e.g. P31946-1) collapsed to canonical accession (P31946).",
        ],
    )


def process_hgnc(path: Path) -> None:
    print("\n=== Processing HGNC complete set ===")
    out_dirs = ensure_dirs(path)

    all_uniprot: set[str] = set()
    all_symbols: set[str] = set()
    all_ensembl: set[str] = set()

    eg_to_sym: set[tuple[str, str]] = set()
    eg_to_up: set[tuple[str, str]] = set()
    sym_to_eg: set[tuple[str, str]] = set()
    sym_to_up: set[tuple[str, str]] = set()
    up_to_eg: set[tuple[str, str]] = set()
    up_to_sym: set[tuple[str, str]] = set()

    n_rows = 0
    with (path / "hgnc_complete_set.txt").open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_rows += 1
            symbol = blank_to_none(row.get("symbol"))
            ensembl = blank_to_none(row.get("ensembl_gene_id"))
            uniprots = split_multi(row.get("uniprot_ids"))

            if symbol:
                all_symbols.add(symbol)
            if ensembl:
                ensembl = strip_ensembl_version(ensembl)
                all_ensembl.add(ensembl)
            for u in uniprots:
                all_uniprot.add(canonical_uniprot(u))

            if symbol and ensembl:
                eg_to_sym.add((ensembl, symbol))
                sym_to_eg.add((symbol, ensembl))
            if symbol and uniprots:
                for u in uniprots:
                    cu = canonical_uniprot(u)
                    sym_to_up.add((symbol, cu))
                    up_to_sym.add((cu, symbol))
            if ensembl and uniprots:
                for u in uniprots:
                    cu = canonical_uniprot(u)
                    eg_to_up.add((ensembl, cu))
                    up_to_eg.add((cu, ensembl))

    print(f"  Read {n_rows:,} HGNC rows")
    print(f"  Unique symbols: {len(all_symbols):,}")
    print(f"  Unique Ensembl genes: {len(all_ensembl):,}")
    print(f"  Unique UniProt IDs: {len(all_uniprot):,}")

    db = "hgnc"
    edge_store = {
        ("ensembl_gene_id", "gene_symbol"): edges_from_pairs(eg_to_sym, "ensembl_gene_id", "gene_symbol", db, db),
        ("ensembl_gene_id", "uniprot_id"): edges_from_pairs(eg_to_up, "ensembl_gene_id", "uniprot_id", db, db),
        ("gene_symbol", "ensembl_gene_id"): edges_from_pairs(sym_to_eg, "gene_symbol", "ensembl_gene_id", db, db),
        ("gene_symbol", "uniprot_id"): edges_from_pairs(sym_to_up, "gene_symbol", "uniprot_id", db, db),
        ("uniprot_id", "ensembl_gene_id"): edges_from_pairs(up_to_eg, "uniprot_id", "ensembl_gene_id", db, db),
        ("uniprot_id", "gene_symbol"): edges_from_pairs(up_to_sym, "uniprot_id", "gene_symbol", db, db),
    }
    id_pools = {
        "ensembl_gene_id": all_ensembl,
        "gene_symbol": all_symbols,
        "uniprot_id": all_uniprot,
    }
    write_outputs(
        out_dirs,
        edge_store,
        id_pools,
        source_label="HGNC",
        extra_notes=[
            "Source file: hgnc_complete_set.txt",
            "Gene symbols are HGNC approved symbols.",
            "Multi-value UniProt fields are split into long-format mappings.",
            "Aliases / previous symbols are ignored for this focused three-ID atlas.",
        ],
    )


def main() -> None:
    uniprot_dir = ARTIFACTS / "UniProt"
    hgnc_dir = ARTIFACTS / "HGNC"
    if not (uniprot_dir / "HUMAN_9606_idmapping.dat.gz").exists():
        raise SystemExit(f"Missing UniProt file under {uniprot_dir}")
    if not (hgnc_dir / "hgnc_complete_set.txt").exists():
        raise SystemExit(f"Missing HGNC file under {hgnc_dir}")

    process_uniprot(uniprot_dir)
    process_hgnc(hgnc_dir)
    print("\n=== Done ===")


if __name__ == "__main__":
    main()
