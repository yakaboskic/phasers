"""
Continuous relative success calculation using ROC-style incremental updates.

Efficiently calculates how RS changes as you select T-I pairs by score, using
incremental count updates rather than full recalculation at each threshold.

Based on ROC-style methodology: sort by score descending, incrementally move
pairs from "not genetically supported" to "genetically supported", track how
RS_cum evolves.

Usage Example:
    >>> from rs_core.continuous import calculate_continuous_rs_incremental
    >>> results = calculate_continuous_rs_incremental(
    ...     pharmaprojects_file='raw/pharmaprojects.tsv',
    ...     indications_file='raw/indications.tsv',
    ...     associations_file='raw/mouse_msigdb.pigean.tsv',
    ...     similarity_matrix_file='raw/minikel.similarity.csv',
    ...     score_column='combined',
    ...     source_name='pigean',
    ...     phase_type='ccat'
    ... )
    >>> # results is a DataFrame with columns: alpha, RS_cum, EA_RS, MLRS
"""

import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Tuple
from itertools import groupby
import logging
import math
from scipy.stats import norm

# Import from existing rs_core modules
from .data_loading import load_and_merge_all_data
from .similarity import calculate_all_similarities
from .phase_mapping import add_phase_numeric_columns, get_unique_phases, validate_phases, PHASE_MAPPING
from .forest_generation import deduplicate_phase_data


def extract_deduplicated_records(
    merged_df: pd.DataFrame,
    phase_type: str,
    score_column: str,
    score_strategy: str = 'max'
) -> pd.DataFrame:
    """
    Extract deduplicated T-I pair records with score and max phase.

    This performs the same deduplication as forest_generation.py, but keeps
    the score column for sorting.

    Args:
        merged_df: Merged pharmaprojects/associations/indications data
        phase_type: One of 'ccat', 'hcat', 'acat'
        score_column: Column name containing association scores
        score_strategy: Strategy for handling scores when deduplicating. One of:
            - 'max': Keep row with highest score
            - 'min': Keep row with lowest score
            - 'avg': Average scores within each ti_uid group

    Returns:
        DataFrame with columns: ti_uid, score, max_phase_num, max_phase_name
        One row per ti_uid, deduplicated by max similarity and max phase
    """
    logging.info(f"Extracting deduplicated records for {phase_type}")

    # # Filter to source if specified
    # if source_name:
    #     # Keep both source-matched rows AND NaN source rows (pharmaprojects-only)
    #     filtered = merged_df[
    #         (merged_df['source'] == source_name) |
    #         (merged_df['source'].isna())
    #     ].copy()
    # else:
    filtered = merged_df.copy()

    # Filter to rows with this phase type
    phase_col = phase_type
    phase_num_col = f'{phase_type}num'

    filtered = filtered[filtered[phase_col].notna()].copy()

    if len(filtered) == 0:
        logging.warning(f"No records found for {phase_type}")
        return pd.DataFrame(columns=['ti_uid', 'score', 'max_phase_num', 'max_phase_name'])

    # Deduplicate by ti_uid using the same logic as forest_generation
    # This ensures consistency with the standard RS calculation
    deduplicated = deduplicate_phase_data(filtered, phase_type, score_column, score_strategy)
    
    # Extract relevant columns
    result = pd.DataFrame({
        'ti_uid': deduplicated['ti_uid'],
        'score': deduplicated[score_column],
        'score_column': score_column,
        'max_phase_num': deduplicated[phase_num_col],
        'max_phase_name': deduplicated[phase_col],
        'target_status': deduplicated['target_status']
    })
    # Filter out records without valid scores or phases
    result = result[
        # result['score'].notna() & # We want to keep records with missing scores because they are from pharmaprojects-only 
        result['max_phase_num'].notna()
    ].copy()

    # Only keep records with actual genetic support status (not meta-statuses)
    result = result[
        ~result['target_status'].isin([
            'indication lacks genetic insight',
            'no indication annotated',
            'no target annotated'
        ])
    ].copy()
    print(result.shape)
    logging.info(f"Extracted {len(result)} deduplicated records")

    return result


def calculate_delta_method_ci_incremental(
    N_G: List[int],
    N_nG: List[int],
    X_G: List[int],
    X_nG: List[int],
    K: int = 4,
    confidence: float = 0.95
) -> Tuple[float, float]:
    """
    Calculate confidence interval for cumulative RS using delta method.

    This is the continuous/incremental version of delta_method_rs_ci from rs_calculation.py.
    Uses the same formula but operates on count arrays instead of forest DataFrames.

    Formula:
    - log_rs_total = sum(log(rs_i)) for phases 0 to K-1
    - var_log_rs_total = sum(var_log_rs_i)
    - CI = exp(log_rs_total ± z * sqrt(var_log_rs_total))
    - var_log_rs_i = (1-p_G)/(N_G*p_G) + (1-p_nG)/(N_nG*p_nG)

    Args:
        N_G: Count arrays with max_phase = i for genetic support group
        N_nG: Count arrays with max_phase = i for non-genetic support group
        X_G: Count arrays advancing beyond phase j for genetic support group
        X_nG: Count arrays advancing beyond phase j for non-genetic support group
        K: Maximum phase (4 for Preclinical=0 to Launched=4)
        confidence: Confidence level (default 0.95 for 95% CI)

    Returns:
        Tuple of (lower_ci, upper_ci)

    Example:
        >>> N_G = [0, 10, 5, 2, 1]
        >>> N_nG = [100, 40, 20, 10, 5]
        >>> X_G = [18, 8, 3, 1]
        >>> X_nG = [75, 50, 15, 5]
        >>> lower, upper = calculate_delta_method_ci_incremental(N_G, N_nG, X_G, X_nG)
        >>> lower < upper
        True
    """
    z = norm.ppf(1 - (1 - confidence) / 2)
    log_rs_total = 0.0
    var_log_rs_total = 0.0

    # Only iterate through phases 0 to K-1 (not phase K = Launched)
    for j in range(0, K):
        # Count that reached AT LEAST phase j
        Ng_total = sum(N_G[k] for k in range(j, K + 1))
        Nng_total = sum(N_nG[k] for k in range(j, K + 1))

        # Count that advanced BEYOND phase j
        Xg = X_G[j]
        Xng = X_nG[j]

        # Check for zero denominators
        if Ng_total == 0 or Nng_total == 0:
            # Can't calculate RS if no targets in either group
            break

        # Calculate success proportions
        p_G = Xg / Ng_total
        p_nG = Xng / Nng_total

        # Check for zero proportions
        if p_nG == 0 or p_G == 0:
            # Can't calculate RS if denominator is zero or log is undefined
            break

        # Calculate relative success
        rs = p_G / p_nG

        if rs <= 0:
            # Can't take log of zero or negative
            break

        # Calculate log(RS) and add to total
        log_rs = np.log(rs)
        log_rs_total += log_rs

        # Calculate variance of log(RS) using delta method
        # Var(log(RS)) = Var(log(p_G)) + Var(log(p_nG))
        # Var(log(p)) ≈ (1-p)/(N*p) for proportion p from N trials
        var_log_rs = (1 - p_G) / (Ng_total * p_G) + (1 - p_nG) / (Nng_total * p_nG)
        var_log_rs_total += var_log_rs

    # Check if we have any variance
    if var_log_rs_total == 0:
        return np.nan, np.nan

    # Calculate confidence interval
    error_margin = z * np.sqrt(var_log_rs_total)
    lower = np.exp(log_rs_total - error_margin)
    upper = np.exp(log_rs_total + error_margin)

    return lower, upper


def build_rs_curve_incremental(
    records_df: pd.DataFrame,
    K: int = 4,
    score_ascending: bool = False
) -> List[Tuple[float, float, float, float, Dict[int, float]]]:
    """
    Build ROC-style RS curve using incremental count updates.

    Algorithm:
    1. Sort records by score (descending by default, or ascending for metrics like p-values)
    2. Initialize: all records in "not genetically supported" (not-G)
    3. For each record (or block of tied scores):
       - Move from not-G to G
       - Update counts: N_G[i]+=1, N_nG[i]-=1, X_G[j]+=1, X_nG[j]-=1
       - Calculate RS_cum = product of RS_j across phases
       - Calculate 95% CI using delta method
       - Record (alpha, RS_cum, lower_ci, upper_ci) where alpha = fraction selected

    Args:
        records_df: DataFrame with columns 'score' and 'max_phase_num'
        K: Maximum phase (4 for Preclinical=0 to Launched=4)
        score_ascending: If True, sort scores ascending (lower is better, e.g., p-values).
                         If False (default), sort descending (higher is better).

    Returns:
        List of (alpha, RS_cum, lower_ci, upper_ci, RS_per_phase, score, N_G, N_nG, X_G, X_nG) tuples where:
        - alpha: fraction of records selected (0 to 1)
        - RS_cum: cumulative relative success
        - lower_ci: lower bound of 95% CI
        - upper_ci: upper bound of 95% CI
        - RS_per_phase: dict mapping phase -> RS at that phase
        - score: score value at this point
        - N_G: list of counts with max_phase=i for genetic support group
        - N_nG: list of counts with max_phase=i for non-genetic support group
        - X_G: list of counts advancing beyond phase j for genetic support group
        - X_nG: list of counts advancing beyond phase j for non-genetic support group

    Note:
        Uses Laplace smoothing to handle zero denominators
    """
    # Sort by score (descending by default, ascending for p-values etc.)
    recs = records_df.sort_values('score', ascending=score_ascending).reset_index(drop=True)

    logging.info(f"Building RS curve for {len(recs)} records, K={K}")

    # Initialize counts
    # N[i] = count with max_phase = i (exact phase, not "at least")
    # X[j] = count with max_phase > j (advanced beyond phase j)
    N_G = [0] * (K + 1)
    N_nG = [0] * (K + 1)
    X_G = [0] * K
    X_nG = [0] * K

    # Start with all records in not-G
    for idx, row in recs.iterrows():
        i = int(row['max_phase_num'])
        N_nG[i] += 1
        for j in range(0, i):  # advanced beyond phases 0..i-1
            X_nG[j] += 1

    logging.info(f"Initial not-G counts: N_nG={N_nG}, X_nG={X_nG}")

    def calculate_rs_cum():
        """Calculate cumulative RS across all phases."""
        RS_cum = 1.0
        rs_per_phase = {}

        for j in range(0, K):  # phases 0 to K-1 (exclude last phase)
            # Count that reached AT LEAST phase j
            # = sum of all records with max_phase >= j
            Ng_total = sum(N_G[k] for k in range(j, K + 1))
            Nng_total = sum(N_nG[k] for k in range(j, K + 1))

            # Count that advanced BEYOND phase j
            Xg = X_G[j]
            Xng = X_nG[j]

            # Laplace smoothing
            pG = (Xg + 0.5) / (Ng_total + 1.0) if Ng_total >= 0 else float('nan')
            pNG = (Xng + 0.5) / (Nng_total + 1.0) if Nng_total >= 0 else float('nan')

            if pNG > 0 and not math.isnan(pG) and not math.isnan(pNG):
                RS_j = pG / pNG
                RS_cum *= RS_j
                rs_per_phase[j] = RS_j
            else:
                rs_per_phase[j] = float('nan')

        return RS_cum, rs_per_phase

    # Build curve by incrementally moving records to G
    selected = 0
    points = []

    # Group by score to handle ties together
    for score, group in recs.groupby('score', sort=False):
        for _, row in group.iterrows():
            i = int(row['max_phase_num'])
            if row['target_status'] == 'genetically supported target': # This means the ti_uid had similarity >= similarity_threshold
                N_G[i] += 1
                N_nG[i] -= 1
                for j in range(0, i):
                    X_G[j] += 1
                    X_nG[j] -= 1
        selected += group.shape[0]
    # for score, group in groupby(zip(recs['score'], recs['max_phase_num'], recs['target_status']), key=lambda x: x[0]):
    #     block = list(group)

    #     # Move entire block from not-G to G
    #     for score_val, max_phase, target_status in block:
    #         i = int(max_phase)
    #         if target_status == 'genetically supported target': # This means the ti_uid had similarity >= similarity_threshold
    #             N_G[i] += 1
    #             N_nG[i] -= 1
    #             for j in range(0, i):
    #                 X_G[j] += 1
    #                 X_nG[j] -= 1

        # selected += len(block)

        # Calculate RS at this point
        RS_cum, rs_per_phase = calculate_rs_cum()
        alpha = selected / len(recs)

        # Calculate confidence intervals using delta method
        lower_ci, upper_ci = calculate_delta_method_ci_incremental(
            N_G, N_nG, X_G, X_nG, K=K, confidence=0.95
        )

        # Store copies of the count lists for debugging
        points.append((alpha, RS_cum, lower_ci, upper_ci, rs_per_phase, score, list(N_G), list(N_nG), list(X_G), list(X_nG)))

    logging.info(f"Generated {len(points)} points on RS curve")

    return points


def calculate_summary_metrics(points: List[Tuple[float, float, float, float, Dict]]) -> Dict[str, float]:
    """
    Calculate scalar summary metrics from RS curve.

    Metrics:
    - EA_RS: Excess Area over RS=1 (average multiplicative lift)
    - MLRS: Mean Log RS (symmetric, robust measure)
    - RS_max: Maximum RS value achieved
    - RS_at_10pct: RS when selecting top 10% of pairs

    Args:
        points: List of (alpha, RS_cum, lower_ci, upper_ci, rs_per_phase) tuples

    Returns:
        Dict with summary metrics
    """
    if len(points) == 0:
        return {
            'EA_RS': float('nan'),
            'MLRS': float('nan'),
            'RS_max': float('nan'),
            'RS_at_10pct': float('nan')
        }

    # Extract alphas and RS values
    alphas = [p[0] for p in points]
    rs_values = [p[1] for p in points]

    # Ensure we have anchor at alpha=0 with RS=1
    if alphas[0] > 0:
        alphas = [0.0] + alphas
        rs_values = [1.0] + rs_values

    # Trapezoidal integration for EA-RS
    ea_rs = sum(
        (alphas[i] - alphas[i-1]) * 0.5 * ((rs_values[i] - 1.0) + (rs_values[i-1] - 1.0))
        for i in range(1, len(alphas))
    )

    # Trapezoidal integration for MLRS
    log_rs_values = [
        math.log(rs) if rs > 0 and not math.isnan(rs) else float('-inf')
        for rs in rs_values
    ]
    mlrs = sum(
        (alphas[i] - alphas[i-1]) * 0.5 * (log_rs_values[i] + log_rs_values[i-1])
        for i in range(1, len(alphas))
        if not math.isinf(log_rs_values[i]) and not math.isinf(log_rs_values[i-1])
    )

    # Max RS
    rs_max = max(rs for rs in rs_values if not math.isnan(rs))

    # RS at 10%
    rs_at_10pct = float('nan')
    for point in points:
        alpha = point[0]
        rs = point[1]
        if alpha >= 0.10:
            rs_at_10pct = rs
            break

    return {
        'EA_RS': ea_rs,
        'MLRS': mlrs,
        'RS_max': rs_max,
        'RS_at_10pct': rs_at_10pct
    }


def calculate_continuous_rs_incremental(
    score_column: str,
    source_name: str,
    phase_type: str = 'ccat',
    similarity_threshold: float = 0.8,
    score_strategy: str = 'max',
    score_ascending: bool = False,
    merged_df: Optional[pd.DataFrame] = None,
    pharmaprojects_file: Optional[str] = None,
    indications_file: Optional[str] = None,
    associations_file: Optional[str] = None,
    similarity_matrix_file: Optional[str] = None,
    similarity_matrix_format: str = 'long',
    no_genetic_insight: bool = True,
    deduplicated_data_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Calculate continuous RS curve using efficient incremental updates.

    This is the main entry point. It:
    1. Loads and merges all data (OR accepts pre-merged data)
    2. Calculates similarities (if not already done)
    3. Deduplicates T-I pairs
    4. Saves deduplicated data (optional)
    5. Builds RS curve using incremental count updates
    6. Calculates summary metrics
    7. Returns points DataFrame

    Args:
        score_column: Column name containing association scores (e.g., 'combined')
        source_name: Source name for filtering (e.g., 'pigean')
        phase_type: Phase type to analyze ('ccat', 'hcat', or 'acat')
        similarity_threshold: Similarity threshold for genetic support (default 0.8)
        score_strategy: Strategy for handling scores when deduplicating. One of:
            - 'max': Keep row with highest score (default)
            - 'min': Keep row with lowest score
            - 'avg': Average scores within each ti_uid group
        score_ascending: If True, sort scores ascending (lower is better, e.g., p-values).
                         If False (default), sort descending (higher is better).
        merged_df: Optional pre-merged DataFrame (if provided, skips loading/merging)
        pharmaprojects_file: Path to pharmaprojects TSV (required if merged_df not provided)
        indications_file: Path to indications TSV (required if merged_df not provided)
        associations_file: Path to associations TSV (required if merged_df not provided)
        similarity_matrix_file: Path to similarity matrix CSV (required if merged_df not provided)
        similarity_matrix_format: Format of similarity matrix ('long' or 'matrix')
        no_genetic_insight: If True, don't filter by genetic insight
        deduplicated_data_file: Optional path to save deduplicated data for inspection

    Returns:
        DataFrame with columns:
        - alpha: Fraction of T-I pairs selected (0 to 1)
        - RS_cum: Cumulative relative success at this selection fraction
        - RS_cum_lower_ci: Lower bound of 95% CI for RS_cum
        - RS_cum_upper_ci: Upper bound of 95% CI for RS_cum
        - RS_phase_0, RS_phase_1, ...: Per-phase RS values
        - EA_RS, MLRS, RS_max, RS_at_10pct: Summary metrics (same in all rows)
    """
    logging.info("=" * 60)
    logging.info("CONTINUOUS RS CALCULATION (INCREMENTAL)")
    logging.info("=" * 60)
    logging.info(f"Score column: {score_column}")
    logging.info(f"Source: {source_name}")
    logging.info(f"Phase type: {phase_type}")
    logging.info(f"Score ascending: {score_ascending} ({'lower is better' if score_ascending else 'higher is better'})")

    # Step 1: Load data (either from merged_df or individual files)
    if merged_df is not None:
        logging.info("Using pre-merged data (skipping load and similarity calculation)")
        logging.info(f"Merged data has {len(merged_df)} rows")

        # Verify required columns exist
        required_cols = ['similarity', 'source', score_column]
        missing_cols = [col for col in required_cols if col not in merged_df.columns]
        if missing_cols:
            raise ValueError(
                f"Merged data is missing required columns: {missing_cols}. "
                f"Available columns: {list(merged_df.columns)}"
            )
    else:
        # Load and merge all data from individual files
        if not all([pharmaprojects_file, indications_file, associations_file, similarity_matrix_file]):
            raise ValueError(
                "Either merged_df OR all of (pharmaprojects_file, indications_file, "
                "associations_file, similarity_matrix_file) must be provided"
            )

        logging.info("Loading and merging data from individual files...")
        merged_df, similarity_matrix_df = load_and_merge_all_data(
            pharmaprojects_file=pharmaprojects_file,
            indications_file=indications_file,
            associations_file=associations_file,
            similarity_matrix_file=similarity_matrix_file,
            similarity_matrix_format=similarity_matrix_format,
            no_genetic_insight=no_genetic_insight,
            association_thresholds=None,  # No filtering - use all data
        )

        logging.info(f"Loaded {len(merged_df)} rows after merge")

        # Step 2: Calculate similarities
        logging.info("Calculating similarities...")
        merged_df = calculate_all_similarities(merged_df, similarity_matrix_df)

    # Step 3: Add phase numeric columns and annotate target status
    merged_df = add_phase_numeric_columns(merged_df)

    # Annotate target_status based on similarity threshold
    merged_df['target_status'] = 'unsupported target'
    merged_df.loc[
        merged_df['similarity'] >= similarity_threshold,
        'target_status'
    ] = 'genetically supported target'

    merged_df.loc[
        merged_df['gene'].isna() | (merged_df['gene'] == ''),
        'target_status'
    ] = 'no target annotated'

    merged_df.loc[
        merged_df['indication_mesh_id'].isna() | (merged_df['indication_mesh_id'] == ''),
        'target_status'
    ] = 'no indication annotated'

    # Step 4: Extract deduplicated records
    dedupped_records = extract_deduplicated_records(
        merged_df,
        phase_type=phase_type,
        score_column=score_column,
        score_strategy=score_strategy
    )

    if len(dedupped_records) == 0:
        raise ValueError(f"No valid deduplicated records found for {source_name}/{phase_type}")

    # Save deduplicated data if requested (sanity check)
    if deduplicated_data_file:
        logging.info(f"Saving deduplicated data to {deduplicated_data_file}")
        dedupped_records.to_csv(deduplicated_data_file, index=False)

    # Step 5: Build RS curve using incremental updates
    logging.info("Building RS curve with incremental updates...")
    K = 4  # Phases 0-4 (Preclinical to Launched)
    points = build_rs_curve_incremental(dedupped_records, K=K, score_ascending=score_ascending)

    # Step 6: Calculate summary metrics
    summary_metrics = calculate_summary_metrics(points)
    logging.info(f"Summary metrics: {summary_metrics}")

    # Step 7: Format output DataFrame
    results_list = []
    for alpha, RS_cum, lower_ci, upper_ci, rs_per_phase, score, N_G, N_nG, X_G, X_nG in points:
        row = {
            'alpha': alpha,
            'RS_cum': RS_cum,
            'RS_cum_lower_ci': lower_ci,
            'RS_cum_upper_ci': upper_ci,
            'source': source_name,
            'phase_type': phase_type,
            'score': score,
        }
        # Add per-phase RS values
        for phase_num, rs_val in rs_per_phase.items():
            row[f'RS_phase_{phase_num}'] = rs_val
        # Add debug counts for each phase (K)
        for k in range(len(N_G)):
            row[f'N_G_{k}'] = N_G[k]
            row[f'N_nG_{k}'] = N_nG[k]
        for k in range(len(X_G)):
            row[f'X_G_{k}'] = X_G[k]
            row[f'X_nG_{k}'] = X_nG[k]
        # Add summary metrics (same for all rows)
        row.update(summary_metrics)
        results_list.append(row)

    results_df = pd.DataFrame(results_list)

    logging.info("=" * 60)
    logging.info("CONTINUOUS RS CALCULATION COMPLETE")
    logging.info(f"Generated {len(results_df)} points on RS curve")
    logging.info(f"EA-RS: {summary_metrics['EA_RS']:.4f}")
    logging.info(f"MLRS: {summary_metrics['MLRS']:.4f}")
    logging.info(f"RS_max: {summary_metrics['RS_max']:.4f}")
    logging.info("=" * 60)

    return results_df


def write_continuous_rs_points(
    results_df: pd.DataFrame,
    output_file: str
) -> None:
    """
    Write continuous RS points to CSV.

    Args:
        results_df: DataFrame from calculate_continuous_rs_incremental()
        output_file: Output CSV file path
    """
    logging.info(f"Writing continuous RS points to {output_file}")
    results_df.to_csv(output_file, index=False)
    logging.info(f"Wrote {len(results_df)} points to {output_file}")
