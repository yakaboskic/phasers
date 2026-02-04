# phasers

**Clinical Phase Relative Success Calculator**

Calculate Relative Success (RS) metrics to evaluate how well genetic evidence predicts drug target success across clinical trial phases.

## Overview

Phasers implements the Relative Success methodology for evaluating genetic support of drug targets. It measures the enrichment of genetically-supported targets among drugs that advance through clinical trial phases (Preclinical → Phase I → Phase II → Phase III → Launched).

**Key Features:**
- Continuous RS curves showing RS vs. selection threshold
- Delta method confidence intervals
- Support for multiple scoring methods (PIGEAN, Open Targets, MAGMA, etc.)
- Interactive Plotly visualizations
- Flexible data loading and filtering

## Installation

Install from GitHub using `uv`:

```bash
uv pip install git+https://github.com/cyakaboski/phasers.git
```

Or using pip:

```bash
pip install git+https://github.com/cyakaboski/phasers.git
```

For development:

```bash
git clone https://github.com/cyakaboski/phasers.git
cd phasers
uv pip install -e .
```

## Quick Start

### Workflow Overview

The typical phasers workflow has three steps:

```
1. MERGE: Combine input files → merged_data.tsv.gz
2. CONTINUOUS: Calculate RS curves → rs_points.csv
3. PLOT: Visualize results → rs_plot.html
```

### Step 1: Merge Data Files

First, merge your input data files and calculate similarities:

```bash
phasers merge \
    --pharmaprojects-file pharmaprojects.tsv \
    --indications-file indications.tsv \
    --associations-file associations.tsv \
    --similarity-matrix-file similarity.csv \
    --output-file merged_data.tsv.gz
```

This creates a reusable merged file that speeds up subsequent analyses.

### Step 2: Calculate Continuous RS

Calculate RS curves for one or more sources and score columns:

```bash
phasers continuous \
    --merged-data-file merged_data.tsv.gz \
    --score-column combined,indirect,direct,overall,drug,pos_p,l2g_share,p_value:ascending \
    --source-name pigean \
    --output-file rs_curve.csv
```

**Key points:**
- **`--score-column`**: List ALL score columns from any source. Columns that don't exist for a source are automatically skipped.
- **`--source-name`**: Can specify multiple sources (comma-separated) to analyze them together.
- Use `column:ascending` suffix for metrics where lower is better (e.g., p-values).

**Real-world example with multiple sources and filters:**

```bash
phasers continuous \
    --merged-data-file merged_data.tsv.gz \
    --score-column combined,indirect,direct,overall,drug,pos_p,l2g_share,p_value:ascending \
    --source-name pigean,falcon,magma,otp \
    --filter "falcon:p_value:<=:1e-8:" \
    --filter "magma:p_value:<=:2.5e-6:" \
    --phase-type ccat \
    --similarity-matrix-threshold 0.8 \
    --score-strategy max \
    --no-genetic-insight \
    --output-file rs_curve.csv
```

This example:
- Analyzes 4 sources: PIGEAN, FALCON, MAGMA, and Open Targets Platform
- Uses all available score columns (each source uses only the columns it has)
- Applies p-value thresholds to FALCON (≤1e-8) and MAGMA (≤2.5e-6)
- Uses the combined clinical category (`ccat`) for phase assessment

### Step 3: Plot Results

Generate interactive HTML plots:

```bash
phasers plot \
    --input-file rs_curve.csv \
    --output-file rs_plot.html \
    --start-alpha 0.001 \
    --max-alpha 0.5
```

## Input Data Formats

### 1. Pharmaprojects File (TSV)

Drug-target-indication clinical trial pipeline data.

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Target gene symbol (e.g., "PCSK9") |
| `indication_mesh_id` | string | MeSH ID for the indication (e.g., "D006937") |
| `ti_uid` | string | Unique target-indication identifier |
| `ccat` | string | Combined clinical phase ("Preclinical", "Phase I", "Phase II", "Phase III", "Launched") |
| `hcat` | string | Historical max phase |
| `acat` | string | Active max phase |
| `succ_p_1` | float | Success indicator: Preclinical → Phase I |
| `succ_1_2` | float | Success indicator: Phase I → Phase II |
| `succ_2_3` | float | Success indicator: Phase II → Phase III |
| `succ_3_a` | float | Success indicator: Phase III → Approval |

**Example:**
```tsv
gene	indication_mesh_id	ti_uid	ccat	succ_p_1	succ_1_2	succ_2_3	succ_3_a
PCSK9	D006937	PCSK9-D006937	Launched	1.0	1.0	1.0	1.0
BRCA1	D001943	BRCA1-D001943	Phase II	1.0	1.0	NaN	NaN
```

### 2. Indications File (TSV)

Indication metadata including genetic insight status.

| Column | Type | Description |
|--------|------|-------------|
| `indication_mesh_id` | string | MeSH ID |
| `genetic_insight` | string | "none" or other value indicating genetic evidence exists |
| `areas` | string | Therapeutic areas (comma-separated) |

**Example:**
```tsv
indication_mesh_id	genetic_insight	areas
D006937	gwas	cardiovascular,metabolic
D001943	gwas	oncology
D003920	none	metabolic
```

### 3. Associations File (TSV)

Gene-indication association scores from various sources.

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Target gene symbol |
| `mesh_id` | string | Association MeSH ID (may be list: `['D001', 'D002']`) |
| `source` | string | Source name (e.g., "pigean", "falcon", "magma") |
| `<score_columns>` | float | Source-specific scores (combined, indirect, p_value, etc.) |

**Example:**
```tsv
gene	mesh_id	source	combined	indirect	direct	p_value
PCSK9	D006937	pigean	0.85	0.72	0.65	NaN
PCSK9	D006937	falcon	NaN	NaN	NaN	1.2e-10
BRCA1	D001943	pigean	0.92	0.88	0.45	NaN
```

### 4. Similarity Matrix File (CSV)

MeSH ID pairwise similarity scores.

**Long format** (default):
```csv
mesh_id_1,mesh_id_2,similarity
D006937,D006937,1.0
D006937,D001943,0.45
D001943,D001943,1.0
```

**Matrix format** (use `--similarity-matrix-format matrix`):
```csv
,D006937,D001943
D006937,1.0,0.45
D001943,0.45,1.0
```

### Common Score Columns by Source

Different association sources provide different score columns. When using `--score-column`, list all columns you want to analyze - columns that don't exist for a source are automatically skipped.

| Source | Score Columns | Sort Order |
|--------|---------------|------------|
| **PIGEAN** | `combined`, `indirect`, `direct` | Descending (higher = better) |
| **Open Targets (otp)** | `overall`, `drug`, `indirect`, `combined` | Descending |
| **FALCON** | `p_value`, `pos_p` | `p_value:ascending`, `pos_p` descending |
| **MAGMA** | `p_value` | Ascending (lower = better) |
| **Open Targets Genetics** | `l2g_share` | Descending |

**Example covering all sources:**
```bash
--score-column combined,indirect,direct,overall,drug,pos_p,l2g_share,p_value:ascending
```

## Python API

```python
from phasers import (
    load_and_merge_data,
    calculate_all_similarities,
    calculate_continuous_rs,
    save_figures
)

# Step 1: Load and merge data
merged_df, similarity_matrix = load_and_merge_data(
    pharmaprojects_file='pharmaprojects.tsv',
    indications_file='indications.tsv',
    associations_file='associations.tsv',
    similarity_matrix_file='similarity.csv'
)

# Step 2: Calculate similarities
merged_df = calculate_all_similarities(merged_df, similarity_matrix)

# Step 3: Calculate RS curve
results = calculate_continuous_rs(
    merged_df=merged_df,
    score_column='combined',
    source_name='pigean',
    phase_type='ccat',
    similarity_threshold=0.8
)

# Results DataFrame contains:
# - alpha: Fraction of T-I pairs selected (0 to 1)
# - RS_cum: Cumulative relative success
# - RS_cum_lower_ci, RS_cum_upper_ci: 95% confidence intervals
# - EA_RS, MLRS, RS_max: Summary metrics
```

## CLI Reference

### `phasers merge`

Merge input data files and calculate similarities.

```
phasers merge [OPTIONS]

Required:
  --pharmaprojects-file FILE    Drug-target-indication clinical trial data (TSV)
  --indications-file FILE       Indication metadata (TSV)
  --associations-file FILE      Gene-indication association scores (TSV)
  --similarity-matrix-file FILE MeSH ID similarity scores (CSV)
  --output-file, -o FILE        Output merged data file (TSV or TSV.gz)

Options:
  --pigean-file FILE            Optional separate PIGEAN associations file
  --similarity-matrix-format    Format: 'long' (default) or 'matrix'
  --no-genetic-insight          Include indications without genetic insight
  --filter FILTER               Apply filter (can repeat)
```

### `phasers continuous`

Calculate continuous RS curve.

```
phasers continuous [OPTIONS]

Required:
  --output-file, -o FILE        Output CSV with RS curve points
  --score-column, -s COLS       Comma-separated score columns from ALL sources
                                (missing columns are skipped per source)
                                Use col:ascending for lower-is-better metrics
  --source-name NAME            Association source(s) to analyze (comma-separated)

Input (choose one):
  --merged-data-file FILE       Pre-merged data file (faster)
  --pharmaprojects-file FILE    + --indications-file + --associations-file + --similarity-matrix-file

Options:
  --phase-type TYPE             Phase type: ccat, hcat, acat (default: ccat)
  --similarity-matrix-threshold Genetic support threshold (default: 0.8)
  --score-strategy STRATEGY     Dedup strategy: max, min, avg (default: avg)
  --filter FILTER               Source-specific filter (can repeat)
                                Format: "source:column:operator:value:"
                                Operators: <, >, <=, >=, ==, <> (outside), >< (inside)
                                Example: --filter "falcon:p_value:<=:1e-8:"
  --no-genetic-insight          Disable genetic insight filter
```

**Filter Examples:**
```bash
--filter "falcon:p_value:<=:1e-8:"      # FALCON p-value ≤ 1e-8
--filter "magma:p_value:<=:2.5e-6:"     # MAGMA p-value ≤ 2.5e-6
--filter "pigean:combined:>=:0.5:"      # PIGEAN combined score ≥ 0.5
```

### `phasers plot`

Plot RS curve from points file.

```
phasers plot [OPTIONS]

Required:
  --input-file, -i FILE         Input CSV with RS points
  --output-file, -o FILE        Output HTML plot

Options:
  --start-alpha FLOAT           Min alpha to plot (default: 0.0)
  --max-alpha FLOAT             Max alpha to plot (default: 1.0)
  --score-column-filter COLS    Filter to specific score columns (comma-separated)
  --source-filter SOURCES       Filter to specific sources (comma-separated)
```

### `phasers summary`

Generate summary statistics.

```
phasers summary [OPTIONS]

Required:
  --input-file, -i FILE         Input CSV with RS points
  --output-file, -o FILE        Output summary CSV
```

## Methodology

The Relative Success methodology implemented in phasers is inspired by the work of Minikel et al., but extended to compute continuous RS curves across all selection thresholds rather than using discrete cutoffs.

> Minikel, E.V., Painter, J.L., Dong, C.C. et al. **Refining the impact of genetic evidence on clinical success.** *Nature* 629, 624–629 (2024). https://doi.org/10.1038/s41586-024-07316-0

### Relative Success (RS)

RS measures the relative probability of clinical advancement for genetically-supported vs. unsupported targets:

```
RS_j = P(advance beyond phase j | genetic support) / P(advance beyond phase j | no genetic support)
```

Cumulative RS across all phases:

```
RS_cum = ∏(RS_j) for j in phases
```

An RS > 1 indicates that genetically-supported targets are more likely to advance through clinical trials.

### Continuous RS Curve

Instead of a single RS value at a fixed threshold, the continuous RS curve shows how RS varies as you select T-I pairs by score:

1. Sort T-I pairs by score (descending for most metrics, ascending for p-values)
2. Incrementally add pairs to the "selected" set
3. Calculate RS at each selection fraction (alpha)

This reveals the optimal threshold and how robust the genetic signal is.

### Summary Metrics

- **EA-RS**: Expected Area under RS curve (integral of RS over alpha)
- **MLRS**: Maximum Likelihood RS (RS at optimal threshold)
- **RS_max**: Maximum RS achieved at any threshold

### Confidence Intervals

Confidence intervals are calculated using the delta method, which propagates uncertainty through the product of per-phase RS values.

## License

MIT License - see LICENSE file for details.

## Citation

If you use phasers in your research, please cite:

**Phasers package:**
```
[Citation for phasers to be added]
```

**Original Relative Success methodology:**
```
Minikel, E.V., Painter, J.L., Dong, C.C. et al. Refining the impact of genetic evidence
on clinical success. Nature 629, 624–629 (2024). https://doi.org/10.1038/s41586-024-07316-0
```

## Contributing

Contributions are welcome! Please open an issue or pull request on GitHub.
