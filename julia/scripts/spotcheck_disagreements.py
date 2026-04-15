#!/usr/bin/env python3
"""
Spot-check disagreements between two similarity sources (ontology vs semantic).

Both merged files share the same (mesh_id, indication_mesh_id) universe — only
the `similarity` column differs. We:
  1. Reduce each file to one row per (mesh_id, indication_mesh_id) -> similarity.
  2. Outer-join the two on that key (so we see pairs only present in one).
  3. Bucket pairs by their (above/below) status at the chosen thresholds:
       - both:          both methods >= their threshold
       - ontology_only: ontology >= ot, semantic <  st (or missing)
       - semantic_only: semantic >= st, ontology <  ot (or missing)
       - neither:       both below threshold (dropped from output to save space)
  4. For each pair, attach indication_mesh_term and a sample of genes that link
     to (mesh_id, indication_mesh_id) so they can be eyeballed.
  5. Write a summary plus top-N disagreement CSVs sorted by absolute similarity gap.

Usage (from any cwd):
    uv run python julia/scripts/spotcheck_disagreements.py \
        --ontology-threshold 0.80 --semantic-threshold 0.80
"""
import argparse
from pathlib import Path

import pandas as pd

JULIA_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = JULIA_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

KEY_COLS = ["mesh_id", "indication_mesh_id"]
LOAD_COLS = ["gene", "mesh_id", "indication_mesh_id", "indication_mesh_term", "similarity"]


def load_pairs(path: str) -> pd.DataFrame:
    print(f"Loading {path} ...", flush=True)
    df = pd.read_csv(
        path, sep="\t",
        usecols=LOAD_COLS,
        dtype={
            "gene": "string",
            "mesh_id": "string",
            "indication_mesh_id": "string",
            "indication_mesh_term": "string",
        },
    )
    df = df.dropna(subset=KEY_COLS + ["similarity"])
    df = df[(df["mesh_id"] != "") & (df["indication_mesh_id"] != "")]
    print(f"  {len(df)} rows after dropping nulls", flush=True)
    return df


def reduce_pairs(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(KEY_COLS, as_index=False)
          .agg(similarity=("similarity", "first"),
               indication_mesh_term=("indication_mesh_term", "first"))
    )


def gene_samples(df: pd.DataFrame, max_genes: int = 5) -> pd.DataFrame:
    return (
        df[KEY_COLS + ["gene"]]
        .dropna()
        .drop_duplicates()
        .groupby(KEY_COLS)["gene"]
        .apply(lambda s: ",".join(sorted(set(s))[:max_genes]))
        .reset_index()
        .rename(columns={"gene": "sample_genes"})
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ontology-file", default=str(JULIA_DIR / "merged_data.tsv.gz"))
    ap.add_argument("--semantic-file", default=str(JULIA_DIR / "merged_cosine_miniLM.tsv.gz"))
    ap.add_argument("--ontology-threshold", type=float, default=0.80)
    ap.add_argument("--semantic-threshold", type=float, default=0.80)
    ap.add_argument("--top-n", type=int, default=200)
    ap.add_argument("--out-prefix", default=None,
                    help="Output prefix (default: results/disagreements_t{ot}_t{st})")
    args = ap.parse_args()

    ot = args.ontology_threshold
    st = args.semantic_threshold
    if args.out_prefix is None:
        prefix = RESULTS_DIR / f"disagreements_t{str(ot).replace('.', '')}_t{str(st).replace('.', '')}"
    else:
        prefix = Path(args.out_prefix)

    ont_df = load_pairs(args.ontology_file)
    sem_df = load_pairs(args.semantic_file)

    print("Reducing to unique (mesh_id, indication_mesh_id) pairs ...", flush=True)
    ont_pairs = reduce_pairs(ont_df).rename(columns={"similarity": "sim_ontology"})
    sem_pairs = reduce_pairs(sem_df).rename(columns={"similarity": "sim_semantic"})

    ont_pairs = ont_pairs.rename(columns={"indication_mesh_term": "indication_mesh_term_ont"})
    sem_pairs = sem_pairs.rename(columns={"indication_mesh_term": "indication_mesh_term_sem"})

    print("Joining ontology vs semantic ...", flush=True)
    merged = ont_pairs.merge(sem_pairs, on=KEY_COLS, how="outer")
    merged["indication_mesh_term"] = (
        merged["indication_mesh_term_ont"].fillna(merged["indication_mesh_term_sem"])
    )
    merged = merged.drop(columns=["indication_mesh_term_ont", "indication_mesh_term_sem"])

    o_above = merged["sim_ontology"].fillna(-1) >= ot
    s_above = merged["sim_semantic"].fillna(-1) >= st
    merged["bucket"] = "neither"
    merged.loc[o_above & s_above, "bucket"] = "both"
    merged.loc[o_above & ~s_above, "bucket"] = "ontology_only"
    merged.loc[~o_above & s_above, "bucket"] = "semantic_only"

    print("Building gene samples ...", flush=True)
    ont_genes = gene_samples(ont_df).rename(columns={"sample_genes": "sample_genes_ontology_side"})
    sem_genes = gene_samples(sem_df).rename(columns={"sample_genes": "sample_genes_semantic_side"})
    merged = merged.merge(ont_genes, on=KEY_COLS, how="left")
    merged = merged.merge(sem_genes, on=KEY_COLS, how="left")

    merged["sim_gap"] = (
        merged["sim_ontology"].fillna(0) - merged["sim_semantic"].fillna(0)
    )

    summary = (
        merged.groupby("bucket")
        .size()
        .rename("n_pairs")
        .reset_index()
        .sort_values("n_pairs", ascending=False)
    )
    summary_path = f"{prefix}_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote {summary_path}")
    print(summary.to_string(index=False))
    print(f"\n(thresholds: ontology >= {ot}, semantic >= {st})")

    cols_out = [
        "mesh_id", "indication_mesh_id", "indication_mesh_term",
        "sim_ontology", "sim_semantic", "sim_gap", "bucket",
        "sample_genes_ontology_side", "sample_genes_semantic_side",
    ]

    def write_top(df, name, by, ascending=False):
        out = df.sort_values(by, ascending=ascending).head(args.top_n)
        path = f"{prefix}_{name}_top{args.top_n}.csv"
        out[cols_out].to_csv(path, index=False)
        print(f"Wrote {path}  ({len(out)} rows)")

    write_top(merged[merged["bucket"] == "ontology_only"], "ontology_only", "sim_ontology")
    write_top(merged[merged["bucket"] == "semantic_only"], "semantic_only", "sim_semantic")
    write_top(merged[merged["bucket"] == "ontology_only"], "ontology_only_largest_gap", "sim_gap")
    write_top(merged[merged["bucket"] == "semantic_only"], "semantic_only_largest_gap", "sim_gap", ascending=True)

    print("\nDone.")


if __name__ == "__main__":
    main()
