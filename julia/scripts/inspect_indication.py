#!/usr/bin/env python3
"""
Inspect what associations have high semantic (and ontology) similarity for a
given indication mesh ID. Used to figure out *what* the embedding model was
actually given as input — the indication mesh term alone? against the
trait_portal label? against the association mesh canonical name?

Usage (from any cwd):
    uv run python julia/scripts/inspect_indication.py D003924
"""
import sys
from pathlib import Path

import pandas as pd

JULIA_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = JULIA_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

INDICATION = sys.argv[1] if len(sys.argv) > 1 else "D003924"  # T2D
TOP_N = 30

COLS = [
    "mesh_id", "indication_mesh_id", "indication_mesh_term",
    "trait_portal", "source", "subsource", "similarity",
]


def load(path):
    print(f"Loading {path} ...", flush=True)
    df = pd.read_csv(
        path, sep="\t",
        usecols=COLS,
        dtype={c: "string" for c in COLS if c != "similarity"},
    )
    df = df[df["indication_mesh_id"] == INDICATION]
    print(f"  {len(df)} rows for indication {INDICATION}", flush=True)
    return df


def top_unique(df, label):
    if df.empty:
        print(f"\n{label}: no rows")
        return None
    term = df["indication_mesh_term"].dropna().iloc[0] if df["indication_mesh_term"].notna().any() else "?"
    print(f"\n=== {label} — indication {INDICATION} ({term}) ===")
    g = (
        df.groupby(["mesh_id", "trait_portal"], dropna=False)
        .agg(similarity=("similarity", "first"),
             source=("source", "first"),
             subsource=("subsource", "first"),
             n_rows=("similarity", "size"))
        .reset_index()
        .sort_values("similarity", ascending=False)
    )
    print(f"\nTop {TOP_N} (mesh_id, trait_portal) by similarity:")
    print(g.head(TOP_N).to_string(index=False))
    print(f"\nBottom {TOP_N} (mesh_id, trait_portal) by similarity:")
    print(g.tail(TOP_N).to_string(index=False))
    return g


def main():
    ont = load(JULIA_DIR / "merged_data.tsv.gz")
    sem = load(JULIA_DIR / "merged_cosine_miniLM.tsv.gz")

    g_ont = top_unique(ont, "ONTOLOGY")
    g_sem = top_unique(sem, "SEMANTIC")

    if g_sem is not None:
        out = RESULTS_DIR / f"inspect_{INDICATION}_semantic.csv"
        g_sem.to_csv(out, index=False)
        print(f"\nWrote {out}")
    if g_ont is not None:
        out = RESULTS_DIR / f"inspect_{INDICATION}_ontology.csv"
        g_ont.to_csv(out, index=False)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
