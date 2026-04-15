# Ontology vs Semantic Similarity — phasers comparison

Comparing two `merged_data` files for phasers continuous-RS:

| file | similarity source | how `similarity` is computed |
|---|---|---|
| `merged_data.tsv.gz` | MeSH ontology distance (existing) | `get_similarity_wrapper` in `src/phasers/core/similarity.py:99` — **max** over the `mesh_id` list |
| `merged_cosine_miniLM.tsv.gz` | sentence-transformers `all-MiniLM-L6-v2` cosine (intern v1) | embeds the **`trait_portal` label** against the **`indication_mesh_term`**, no max-over-list aggregation |

Both files share the same row universe and the same `mesh_id` lists; only the
`similarity` column differs.

## Layout

```
julia/
├── README.md
├── merged_data.tsv.gz             # ontology similarity
├── merged_cosine_miniLM.tsv.gz    # semantic v1 (intern)
├── scripts/
│   ├── run_phasers.sh             # interactive single-file runner
│   ├── sweep_ea_rs.sh             # ontology@0.8 + semantic sweep
│   ├── parse_ea_rs.py             # parses sweep logs -> ea_rs_sweep.csv
│   ├── analyze_gst_overlap.py     # GST count + Jaccard sweep
│   ├── spotcheck_disagreements.py # per-pair disagreements at chosen thresholds
│   └── inspect_indication.py      # per-indication top-N association inspector
├── results/                       # all generated CSVs and logs
└── archive/                       # earlier intern runs
```

All scripts resolve paths from their own location, so they can be invoked from
any working directory:

```bash
uv run python julia/scripts/analyze_gst_overlap.py
uv run python julia/scripts/spotcheck_disagreements.py --ontology-threshold 0.8 --semantic-threshold 0.8
uv run python julia/scripts/inspect_indication.py D003924
julia/scripts/sweep_ea_rs.sh
uv run python julia/scripts/parse_ea_rs.py
```

## What the analysis found (semantic v1)

Anchor: **ontology @ 0.8 → EA-RS = 0.8031** (2,176 unique GST genes).

1. **No semantic threshold matches the ontology anchor on EA-RS.** Best
   semantic candidates plateau at EA-RS ≈ 0.65–0.66 (thresholds 0.30 and 0.70).
   Matched-count threshold (semantic @ 0.55) gives EA-RS = 0.5607. See
   `results/ea_rs_sweep.csv`.
2. **GST overlap at matched count is high on genes (Jaccard 0.90), low on
   gene–indication pairs (Jaccard 0.68).** The two methods generally pick the
   same target *genes* but disagree on *which indications* those genes are
   supported for. See `results/gst_overlap_summary.csv`.
3. **The ontology gets free wins from direct mesh-set membership.** Many
   ontology pairs score 1.0 because the indication's mesh ID is literally an
   element of the association's `mesh_id` list (e.g. Type 2 Diabetes inside
   `['D003920','D003922','D003924','D006946',...]`). Semantic v1 has no
   information about set membership and scores those same pairs at 0.05–0.30.
   See `results/disagreements_t08_t08_ontology_only_top200.csv`.
4. **Semantic recovers a small slice that ontology has no opinion on**, where
   the ontology returns sentinel `-4.0` (no path). E.g. Huntington Disease vs
   `D020271`. Only 39 such pairs at 0.8/0.8.
5. **Semantic v1 was built against `trait_portal`, not the mesh canonical
   name.** Diagnostic: the pair `mesh_id=['D003924']`, `trait_portal='T2D'`
   gets cosine 0.865 against indication "Diabetes Mellitus, Type 2". If the
   embedding had used the mesh canonical of `D003924` ("Diabetes Mellitus,
   Type 2"), it would have been an exact-string self-match → cosine ≈ 1.0.
   The ranking in `results/inspect_D003924_semantic.csv` is dominated by
   `trait_portal` token overlap (`T2D`, `T2D_dom`, `T2D_age-related`,
   `gcat_trait_type_2_diabetes_mellitus`), not by structured mesh terms.

## Next steps for the semantic embedding rebuild

The fix that lets us cleanly attribute the remaining gap to "embedding model"
vs "input text quality":

### 1. Rebuild `merged_cosine_*.tsv.gz` with max-over-mesh-list aggregation

For each row, the similarity should be:

```
similarity = max(
    cosine(embed(indication_mesh_term), embed(canonical_term(m)))
    for m in mesh_id_list
)
```

This mirrors `get_similarity_wrapper` in `src/phasers/core/similarity.py:99`
(see lines 130–147 — the ontology version iterates the same list and takes
the max). The new file should:

- Use **mesh canonical names** on the association side, not `trait_portal`.
- Aggregate via **max over the `mesh_id` list** (same as the ontology code).
- Preserve the same row count and schema as `merged_data.tsv.gz` so all the
  scripts in this folder run unchanged against the new file.

Once that file exists (suggested name: `merged_cosine_miniLM_v2.tsv.gz`),
re-run the same sweep against it. If the EA-RS gap closes, the remaining
difference between v2 and ontology is purely the embedding distance vs the
ontology distance — clean apples-to-apples.

### 2. (Optional) Try a biomedical embedding model

`all-MiniLM-L6-v2` is general-purpose and doesn't know that `T2D` ≈ "Diabetes
Mellitus, Type 2" or that `CKD` ≈ "Chronic Kidney Disease". Worth retrying
the v2 build with a biomedical-tuned encoder:

- `cambridgeltl/SapBERT-from-PubMedBERT-fulltext`
- `pritamdeka/S-PubMedBert-MS-MARCO`
- `NeuML/pubmedbert-base-embeddings`

### 3. (Optional) A v3 that *also* uses `trait_portal` as auxiliary context

If after v2 there's still a gap and the disagreements look like the model is
missing study-specific context (e.g. an ontology hit on a parent term that a
canonical-name embedding misses), try a v3 that takes the max over both the
mesh canonical names *and* the `trait_portal` text. The intern's v1 used only
`trait_portal`; v2 uses only mesh canonical names; v3 would use both.

### Suggested validation workflow once each new file exists

```bash
# 1. EA-RS sweep against the new file (edit SEMANTIC_FILE in sweep_ea_rs.sh)
julia/scripts/sweep_ea_rs.sh
uv run python julia/scripts/parse_ea_rs.py

# 2. GST count + Jaccard
uv run python julia/scripts/analyze_gst_overlap.py \
    --semantic-file julia/merged_cosine_miniLM_v2.tsv.gz

# 3. Disagreement spot-check at matched thresholds
uv run python julia/scripts/spotcheck_disagreements.py \
    --semantic-file julia/merged_cosine_miniLM_v2.tsv.gz \
    --ontology-threshold 0.8 --semantic-threshold 0.8

# 4. Per-indication sanity checks
uv run python julia/scripts/inspect_indication.py D003924   # T2D
uv run python julia/scripts/inspect_indication.py D006816   # Huntington
uv run python julia/scripts/inspect_indication.py D008175   # Lung Cancer
```

## Reference: ontology similarity code

`src/phasers/core/similarity.py`:

- `get_similarity` (line 30) — single-pair lookup against `similarity_matrix`,
  returns sentinel codes `-1` (both null), `-2` (no indication), `-3` (no
  association), `-4` (mesh ID not in matrix).
- `get_similarity_wrapper` (line 99) — handles `mesh_id` as a list (or a
  string-encoded list) and **takes the max** across all elements.
- `calculate_all_similarities` (line 154) — applies the wrapper row-wise
  across the merged dataframe.

Confirmed by reading the source: ontology `merged_data.tsv.gz` was built with
exactly this max-over-list aggregation.
