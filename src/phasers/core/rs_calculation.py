"""
Core relative success calculation with confidence intervals.

Implements:
- Cumulative count calculations
- P(G) and P(S) proportion calculations
- Relative success ratios
- Wilson confidence intervals for proportions
- Delta method confidence intervals for overall RS

Reference: relative_success.py lines 97-137, 273-408
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.proportion import proportion_confint
from typing import Tuple, Dict


def calculate_rs_for_forest(forest_df: pd.DataFrame, phase_type: str) -> pd.DataFrame:
    """
    Calculate relative success metrics for forest data.

    Adds columns:
    - total_targets_cumsum: Reverse cumsum from Launched to Preclinical
    - number_targets_with_genetic_support_cumsum: Reverse cumsum
    - number_targets_without_genetic_support_cumsum: Difference
    - P(G): Proportion with genetic support
    - P(S), P(S)_g, P(S)_ng: Proportion advancing to next phase
    - relative_success: Ratio of P(S)_g / P(S)_ng
    - Various confidence interval columns

    Args:
        forest_df: Forest data for one phase type
        phase_type: One of 'ccat', 'hcat', 'acat'

    Returns:
        Forest DataFrame with added RS calculation columns

    Reference: relative_success.py:273-384

    Note:
        Handles many edge cases:
        - Missing next phase (sets values to NaN)
        - Zero denominators (sets RS to NaN)
        - Zero success proportions (sets RS to NaN)

        Cumulative sums are reverse (Launched to Preclinical) because we want
        to know "of all targets that reached at least this phase, how many
        had genetic support?"

    Example:
        >>> forest = pd.DataFrame({
        ...     'ccatnum': [0, 1, 2, 3, 4],
        ...     'ccat': ['Preclinical', 'Phase I', 'Phase II', 'Phase III', 'Launched'],
        ...     'total_targets': [100, 50, 20, 10, 5],
        ...     'number_targets_with_genetic_support': [60, 35, 15, 8, 4]
        ... })
        >>> result = calculate_rs_for_forest(forest, 'ccat')
        >>> 'relative_success' in result.columns
        True
    """
    # Sort by phase number
    forest_df = forest_df.sort_values(by=[f'{phase_type}num'])

    # Calculate reverse cumulative sums (from Launched back to Preclinical)
    # This gives us "of all targets that reached at least this phase..."
    forest_df['total_targets_cumsum'] = forest_df.iloc[::-1]['total_targets'].cumsum()[::-1]
    forest_df['number_targets_with_genetic_support_cumsum'] = \
        forest_df.iloc[::-1]['number_targets_with_genetic_support'].cumsum()[::-1]
    forest_df['number_targets_without_genetic_support_cumsum'] = (
        forest_df['total_targets_cumsum'] - forest_df['number_targets_with_genetic_support_cumsum']
    )

    # Calculate P(G): Proportion of targets with genetic support
    try:
        forest_df['P(G)'] = (
            forest_df['number_targets_with_genetic_support_cumsum'] /
            forest_df['total_targets_cumsum']
        )
        forest_df['lower_ci'], forest_df['upper_ci'] = proportion_confint(
            forest_df['number_targets_with_genetic_support_cumsum'],
            forest_df['total_targets_cumsum'],
            alpha=0.05,
            method='wilson'
        )
    except ZeroDivisionError:
        print(f"Warning: Zero division error for {phase_type} forest data in P(G) calculation")
        forest_df['P(G)'] = np.nan
        forest_df['lower_ci'] = np.nan
        forest_df['upper_ci'] = np.nan

    # Calculate P(S) and relative success for each phase
    for i, row in forest_df.iterrows():
        phase_label = f'{phase_type}num'
        phase_num = row[phase_label]
        next_phase_num = phase_num + 1

        # Handle Launched phase (no next phase)
        if phase_num == 4:
            forest_df.loc[i, 'P(S)'] = np.nan
            forest_df.loc[i, 'P(S)_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng'] = np.nan
            forest_df.loc[i, 'P(S)_ng_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g'] = np.nan
            forest_df.loc[i, 'P(S)_g_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g_upper_ci'] = np.nan
            forest_df.loc[i, 'relative_success'] = np.nan
            forest_df.loc[i, 'relative_success_lower_ci'] = np.nan
            forest_df.loc[i, 'relative_success_upper_ci'] = np.nan
            continue

        # Handle missing next phase
        elif len(forest_df.loc[forest_df[phase_label] == next_phase_num]) == 0:
            forest_df.loc[i, 'P(S)'] = np.nan
            forest_df.loc[i, 'P(S)_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng'] = np.nan
            forest_df.loc[i, 'P(S)_ng_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g'] = np.nan
            forest_df.loc[i, 'P(S)_g_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g_upper_ci'] = np.nan
            forest_df.loc[i, 'relative_success'] = np.nan
            forest_df.loc[i, 'relative_success_lower_ci'] = np.nan
            forest_df.loc[i, 'relative_success_upper_ci'] = np.nan
            continue

        # Handle zero total targets
        elif row['total_targets_cumsum'] == 0:
            forest_df.loc[i, 'P(S)'] = np.nan
            forest_df.loc[i, 'P(S)_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng'] = np.nan
            forest_df.loc[i, 'P(S)_ng_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_ng_upper_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g'] = np.nan
            forest_df.loc[i, 'P(S)_g_lower_ci'] = np.nan
            forest_df.loc[i, 'P(S)_g_upper_ci'] = np.nan
            forest_df.loc[i, 'relative_success'] = np.nan
            forest_df.loc[i, 'relative_success_lower_ci'] = np.nan
            forest_df.loc[i, 'relative_success_upper_ci'] = np.nan
            continue

        # Calculate P(S): Overall proportion advancing to next phase
        else:
            next_phase_total = forest_df.loc[
                forest_df[phase_label] == next_phase_num,
                'total_targets_cumsum'
            ].values[0]

            forest_df.loc[i, 'P(S)'] = next_phase_total / row['total_targets_cumsum']
            forest_df.loc[i, 'P(S)_lower_ci'], forest_df.loc[i, 'P(S)_upper_ci'] = \
                proportion_confint(
                    next_phase_total,
                    row['total_targets_cumsum'],
                    alpha=0.05,
                    method='wilson'
                )

            # Calculate P(S)_ng: Proportion without genetic support advancing
            next_phase_without = forest_df.loc[
                forest_df[phase_label] == next_phase_num,
                'number_targets_without_genetic_support_cumsum'
            ].values[0]

            forest_df.loc[i, 'P(S)_ng'] = \
                next_phase_without / row['number_targets_without_genetic_support_cumsum']
            forest_df.loc[i, 'P(S)_ng_lower_ci'], forest_df.loc[i, 'P(S)_ng_upper_ci'] = \
                proportion_confint(
                    next_phase_without,
                    row['number_targets_without_genetic_support_cumsum'],
                    alpha=0.05,
                    method='wilson'
                )

            # Handle zero genetic support counts
            if (row['number_targets_with_genetic_support_cumsum'] == 0 or
                row['number_targets_without_genetic_support_cumsum'] == 0):
                forest_df.loc[i, 'P(S)_g'] = np.nan
                forest_df.loc[i, 'P(S)_g_lower_ci'] = np.nan
                forest_df.loc[i, 'P(S)_g_upper_ci'] = np.nan
                forest_df.loc[i, 'relative_success'] = np.nan
                forest_df.loc[i, 'relative_success_lower_ci'] = np.nan
                forest_df.loc[i, 'relative_success_upper_ci'] = np.nan
                continue

            # Calculate P(S)_g: Proportion with genetic support advancing
            next_phase_with = forest_df.loc[
                forest_df[phase_label] == next_phase_num,
                'number_targets_with_genetic_support_cumsum'
            ].values[0]

            forest_df.loc[i, 'P(S)_g'] = \
                next_phase_with / row['number_targets_with_genetic_support_cumsum']
            forest_df.loc[i, 'P(S)_g_lower_ci'], forest_df.loc[i, 'P(S)_g_upper_ci'] = \
                proportion_confint(
                    next_phase_with,
                    row['number_targets_with_genetic_support_cumsum'],
                    alpha=0.05,
                    method='wilson'
                )

            # Calculate relative success ratios
            ratio1 = next_phase_with / row['number_targets_with_genetic_support_cumsum']
            ratio1_lower, ratio1_upper = proportion_confint(
                next_phase_with,
                row['number_targets_with_genetic_support_cumsum'],
                alpha=0.05,
                method='wilson'
            )

            ratio2 = next_phase_without / row['number_targets_without_genetic_support_cumsum']
            ratio2_lower, ratio2_upper = proportion_confint(
                next_phase_without,
                row['number_targets_without_genetic_support_cumsum'],
                alpha=0.05,
                method='wilson'
            )

            # Calculate relative success: RS = P(S)_g / P(S)_ng
            if ratio2 == 0:
                forest_df.loc[i, 'relative_success'] = np.nan
            else:
                forest_df.loc[i, 'relative_success'] = ratio1 / ratio2

            # Calculate confidence intervals for RS
            if ratio2_lower == 0:
                forest_df.loc[i, 'relative_success_upper_ci'] = np.nan
            else:
                forest_df.loc[i, 'relative_success_upper_ci'] = ratio1_upper / ratio2_lower

            if ratio2_upper == 0:
                forest_df.loc[i, 'relative_success_lower_ci'] = np.nan
            else:
                forest_df.loc[i, 'relative_success_lower_ci'] = ratio1_lower / ratio2_upper

    # Reset index
    forest_df = forest_df.reset_index(drop=True)

    return forest_df


def delta_method_rs_ci(forest_df: pd.DataFrame, confidence: float = 0.95) -> Tuple[float, float]:
    """
    Calculate confidence interval for overall RS using delta method.

    Calculates log(RS) and variance across all phases, then exponentiates
    to get the confidence interval for the product of RS values.

    Formula:
    - log_rs_total = sum(log(rs_i)) for phases 0-3
    - var_log_rs_total = sum(var_log_rs_i)
    - CI = exp(log_rs_total ± z * sqrt(var_log_rs_total))

    Args:
        forest_df: Forest DataFrame with RS calculations
        confidence: Confidence level (default 0.95 for 95% CI)

    Returns:
        Tuple of (lower_ci, upper_ci)

    Reference: relative_success.py:97-137

    Note:
        Only uses phases 0-3 (Preclinical to Phase III), not Launched (phase 4)
        because there's no phase after Launched to calculate advancement.

        Breaks early if:
        - Any denominator is zero (can't calculate RS)
        - Any RS is zero (log undefined)
        - Any P(S)_ng is zero (can't calculate variance)

    Example:
        >>> forest = pd.DataFrame({
        ...     'relative_success': [1.5, 1.3, 1.2, 1.1, np.nan],
        ...     'number_targets_with_genetic_support_cumsum': [60, 35, 15, 8, 4],
        ...     'number_targets_without_genetic_support_cumsum': [40, 15, 5, 2, 1]
        ... })
        >>> forest.loc[0, 'number_targets_with_genetic_support_cumsum'] = 60
        >>> lower, upper = delta_method_rs_ci(forest)
        >>> lower < upper
        True
    """
    z = norm.ppf(1 - (1 - confidence) / 2)
    log_rs_total = 0
    var_log_rs_total = 0

    # Only iterate through phases 0-3 (not phase 4 = Launched)
    for row_idx in range(5):
        if row_idx == 4:
            # Skip Launched phase (no next phase to advance to)
            continue

        # Get counts
        N_G = forest_df.iloc[row_idx]['number_targets_with_genetic_support_cumsum']
        N_noG = forest_df.iloc[row_idx]['number_targets_without_genetic_support_cumsum']

        # Get number advancing to next phase
        X_G = forest_df.iloc[row_idx + 1]['number_targets_with_genetic_support_cumsum']
        X_noG = forest_df.iloc[row_idx + 1]['number_targets_without_genetic_support_cumsum']

        # Check for zero denominators
        if N_G == 0 or N_noG == 0:
            # Can't calculate RS if no targets in either group
            break

        # Calculate success proportions
        p_G = X_G / N_G  # Proportion with genetic support advancing
        p_noG = X_noG / N_noG  # Proportion without genetic support advancing

        # Check for zero proportions
        if p_noG == 0:
            # Can't calculate RS if denominator is zero
            break

        # Calculate relative success
        rs = p_G / p_noG

        if rs == 0:
            # Can't take log of zero
            break

        # Calculate log(RS) and add to total
        log_rs = np.log(rs)
        log_rs_total += log_rs

        # Calculate variance of log(RS) using delta method
        # Var(log(RS)) = Var(log(p_G)) + Var(log(p_noG))
        # Var(log(p)) ≈ (1-p)/(N*p) for proportion p from N trials
        var_log_rs = (1 - p_G) / (N_G * p_G) + (1 - p_noG) / (N_noG * p_noG)
        var_log_rs_total += var_log_rs

    # Check if we have any variance
    if var_log_rs_total == 0:
        return np.nan, np.nan

    # Calculate confidence interval
    error_margin = z * np.sqrt(var_log_rs_total)
    lower = np.exp(log_rs_total - error_margin)
    upper = np.exp(log_rs_total + error_margin)

    return lower, upper


def calculate_overall_rs(forest_dict: Dict[str, pd.DataFrame]) -> Dict[Tuple[str, str], float]:
    """
    Calculate overall RS for all phase types with confidence intervals.

    For each phase_type (ccat, hcat, acat):
    - Calculate product of relative_success for phases 0-3
    - Calculate delta method confidence interval

    Args:
        forest_dict: Dict mapping phase_type to forest DataFrame

    Returns:
        Dict with keys (phase_type, metric_type) and values:
            - (phase_type, 'exact'): Point estimate
            - (phase_type, 'lower_ci'): Lower confidence bound
            - (phase_type, 'upper_ci'): Upper confidence bound

    Reference: relative_success.py:388-397

    Example:
        >>> forests = {
        ...     'ccat': pd.DataFrame({
        ...         'relative_success': [1.5, 1.3, 1.2, 1.1, np.nan],
        ...         'number_targets_with_genetic_support_cumsum': [60, 35, 15, 8, 4],
        ...         'number_targets_without_genetic_support_cumsum': [40, 15, 5, 2, 1]
        ...     })
        ... }
        >>> rs_results = calculate_overall_rs(forests)
        >>> ('ccat', 'exact') in rs_results
        True
    """
    relative_success_all = {}

    for phase_type, forest in forest_dict.items():
        # Calculate product of RS for phases 0-3 (skipna=False means NaN propagates)
        # Only include first 4 phases (0-3), not Launched (4)
        relative_success_all[(phase_type, 'exact')] = \
            forest['relative_success'].iloc[:4].product(skipna=False)

        # Calculate confidence interval using delta method
        lower, upper = delta_method_rs_ci(forest, confidence=0.95)
        relative_success_all[(phase_type, 'lower_ci')] = lower
        relative_success_all[(phase_type, 'upper_ci')] = upper

    return relative_success_all
