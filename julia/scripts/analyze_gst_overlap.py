#!/usr/bin/env python3
"""
Compare genetically-supported-target (GST) sets between two merged-data files
that use different similarity sources (e.g., ontology vs. semantic embedding).

For each file, sweep similarity thresholds and count:
  - rows above threshold
  - unique GST genes
  - unique (gene, indication_mesh_id) pairs

Then anchor on a chosen ontology threshold (default 0.8), find the semantic
threshold whose unique-gene count is closest, and report Jaccard similarity
between the two GST sets (both genes and gene-indication pairs).

Usage (from any cwd):
    uv run python julia/scripts/analyze_gst_overlap.py
"""
import argparse
from pathlib import Path

import pandas as pd

JULIA_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = JULIA_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

DEFAULT_THRESHOLDS = [
    0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
    0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
]


def load_min(path):
    print(f"Loading {path} ...", flush=True)
    df = pd.read_csv(
        path, sep="\t",
        usecols=["gene", "indication_mesh_id", "similarity"],
        dtype={"gene": "string", "indication_mesh_id": "string"},
    )
    n0 = len(df)
    df = df.dropna(subset=["gene", "indication_mesh_id", "similarity"])
    df = df[(df["gene"] != "") & (df["indication_mesh_id"] != "")]
    print(f"  {n0} -> {len(df)} rows after dropping null gene/indication/similarity",
          flush=True)
    return df.reset_index(drop=True)


def gst_sets(df, threshold):
    sub = df[df["similarity"] >= threshold]
    genes = set(sub["gene"].unique())
    pairs = set(zip(sub["gene"], sub["indication_mesh_id"]))
    return genes, pairs, len(sub)


def jaccard(a, b):
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ontology-file", default=str(JULIA_DIR / "merged_data.tsv.gz"))
    ap.add_argument("--semantic-file", default=str(JULIA_DIR / "merged_cosine_miniLM.tsv.gz"))
    ap.add_argument("--anchor-threshold", type=float, default=0.8)
    ap.add_argument("--out-sweep", default=str(RESULTS_DIR / "gst_count_sweep.csv"))
    ap.add_argument("--out-summary", default=str(RESULTS_DIR / "gst_overlap_summary.csv"))
    args = ap.parse_args()

    ont = load_min(args.ontology_file)
    sem = load_min(args.semantic_file)

    rows = []
    for t in DEFAULT_THRESHOLDS:
        for label, df in [("ontology", ont), ("semantic", sem)]:
            g, p, n = gst_sets(df, t)
            rows.append({
                "source": label,
                "threshold": t,
                "n_rows": n,
                "n_unique_genes": len(g),
                "n_unique_gene_indication_pairs": len(p),
            })
    sweep = pd.DataFrame(rows)
    sweep.to_csv(args.out_sweep, index=False)
    print(f"\nWrote {args.out_sweep}")
    print(sweep.to_string(index=False))

    ont_genes, ont_pairs, _ = gst_sets(ont, args.anchor_threshold)
    target_n = len(ont_genes)
    print(f"\nAnchor: ontology @ {args.anchor_threshold} -> {target_n} unique GST genes")

    sem_sweep = sweep[sweep["source"] == "semantic"].copy()
    sem_sweep["diff"] = (sem_sweep["n_unique_genes"] - target_n).abs()
    best = sem_sweep.sort_values("diff").iloc[0]
    matched_t = float(best["threshold"])
    print(f"Closest-count semantic threshold: {matched_t} "
          f"-> {int(best['n_unique_genes'])} unique GST genes "
          f"(|diff|={int(best['diff'])})")

    sem_genes, sem_pairs, _ = gst_sets(sem, matched_t)
    sem_genes_same, sem_pairs_same, _ = gst_sets(sem, args.anchor_threshold)

    def row(comparison, og, sg, op, sp):
        return {
            "comparison": comparison,
            "jaccard_genes": jaccard(og, sg),
            "jaccard_pairs": jaccard(op, sp),
            "n_ont_genes": len(og),
            "n_sem_genes": len(sg),
            "n_intersect_genes": len(og & sg),
            "n_union_genes": len(og | sg),
            "n_ont_pairs": len(op),
            "n_sem_pairs": len(sp),
            "n_intersect_pairs": len(op & sp),
            "n_union_pairs": len(op | sp),
        }

    summary = pd.DataFrame([
        row(f"ontology@{args.anchor_threshold} vs semantic@{matched_t} (matched count)",
            ont_genes, sem_genes, ont_pairs, sem_pairs),
        row(f"ontology@{args.anchor_threshold} vs semantic@{args.anchor_threshold} (same threshold)",
            ont_genes, sem_genes_same, ont_pairs, sem_pairs_same),
    ])
    summary.to_csv(args.out_summary, index=False)
    print(f"\nWrote {args.out_summary}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
