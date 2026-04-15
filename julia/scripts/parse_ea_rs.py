#!/usr/bin/env python3
"""
Parse phasers continuous logs from sweep_ea_rs.sh and build a comparison table.

Reports the semantic threshold whose EA-RS is closest to the ontology anchor
(ontology @ 0.8 by default).

Usage (from any cwd):
    uv run python julia/scripts/parse_ea_rs.py
"""
import re
import sys
from pathlib import Path

import pandas as pd

JULIA_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = JULIA_DIR / "results"
DEFAULT_LOG_DIR = RESULTS_DIR / "sweep_logs"

METRIC_KEYS = [
    "EA-RS", "EA-RS-self", "RS_full", "MLRS",
    "RS_max", "RS_at_10pct",
    "RS_max_after_1pct", "RS_max_after_5pct", "RS_max_after_10pct",
]

LOG_RE = re.compile(r"rs_run_(ontology|semantic)_(\d+p\d+)\.log")
GST_RE = re.compile(
    r"Genetically supported targets at similarity_threshold=([\d.]+):\s*"
    r"(\d+) rows,\s*(\d+) unique genes"
)


def parse_log(path: Path) -> dict:
    m = LOG_RE.search(path.name)
    if not m:
        return {}
    label = m.group(1)
    threshold = float(m.group(2).replace("p", "."))
    text = path.read_text()

    out = {"source": label, "threshold": threshold, "log_file": path.name}
    for key in METRIC_KEYS:
        mm = re.search(rf"{re.escape(key)}:\s*(-?[\d.]+)", text)
        if mm:
            out[key.replace("-", "_")] = float(mm.group(1))
    g = GST_RE.search(text)
    if g:
        out["n_gst_rows"] = int(g.group(2))
        out["n_gst_genes"] = int(g.group(3))
    return out


def main():
    log_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_LOG_DIR
    rows = [parse_log(p) for p in sorted(log_dir.glob("rs_run_*.log"))]
    rows = [r for r in rows if r]
    if not rows:
        print(f"No matching logs in {log_dir}", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows).sort_values(["source", "threshold"]).reset_index(drop=True)
    out = RESULTS_DIR / "ea_rs_sweep.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}\n")
    print(df.to_string(index=False))

    anchor = df[(df["source"] == "ontology") & (df["threshold"] == 0.80)]
    if anchor.empty:
        print("\nNo ontology @ 0.80 row found; cannot anchor.")
        return
    anchor_ea = float(anchor["EA_RS"].iloc[0])
    print(f"\nAnchor: ontology @ 0.80 -> EA_RS = {anchor_ea:.4f}")

    sem = df[df["source"] == "semantic"].copy()
    sem["ea_rs_diff"] = (sem["EA_RS"] - anchor_ea).abs()
    sem_sorted = sem.sort_values("ea_rs_diff")
    best = sem_sorted.iloc[0]
    print(
        f"Closest semantic threshold by EA_RS: {best['threshold']:.2f} "
        f"(EA_RS={best['EA_RS']:.4f}, |diff|={best['ea_rs_diff']:.4f})"
    )

    print("\nSemantic candidates ranked by |EA_RS - anchor|:")
    cols = ["threshold", "EA_RS", "ea_rs_diff", "RS_full", "MLRS", "n_gst_genes"]
    cols = [c for c in cols if c in sem_sorted.columns]
    print(sem_sorted[cols].to_string(index=False))


if __name__ == "__main__":
    main()
