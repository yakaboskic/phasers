"""
Forest data generation for relative success analysis.

Handles:
- Deduplication of ti_uid pairs by max similarity and phase
- Grouping by phase to count targets with/without genetic support
- Missing phase handling

Reference: relative_success.py lines 87-95, 213-270
"""

import pandas as pd
from typing import Set, Optional
from .phase_mapping import PHASE_MAPPING_REVERSE


def get_forest_data(phase_type: str, phase_data: pd.DataFrame) -> pd.DataFrame:
    """
    Generate forest data for a single phase.

    Groups data by phase and calculates:
    - total_targets: Count of targets (excluding certain statuses)
    - number_targets_with_genetic_support: Count where target_status == 'genetically supported target'

    Args:
        phase_type: One of 'ccat', 'hcat', 'acat'
        phase_data: DataFrame for a specific phase, already deduplicated

    Returns:
        DataFrame with phase statistics

    Reference: relative_success.py:87-95

    Note:
        Excludes target_status values:
        - 'indication lacks genetic insight'
        - 'no indication annotated'
        - 'no target annotated'

        These are excluded because they don't represent actual drug targets
        that could have genetic support.

    Example:
        >>> data = pd.DataFrame({
        ...     'ccatnum': [0, 0, 0],
        ...     'ccat': ['Preclinical', 'Preclinical', 'Preclinical'],
        ...     'target_status': ['genetically supported target',
        ...                       'unsupported target',
        ...                       'genetically supported target']
        ... })
        >>> forest = get_forest_data('ccat', data)
        >>> forest['total_targets'].iloc[0]
        3
        >>> forest['number_targets_with_genetic_support'].iloc[0]
        2
    """
    # Filter out targets with certain statuses
    # These don't represent actual evaluable targets
    forest = phase_data[
        ~phase_data['target_status'].isin([
            'indication lacks genetic insight',
            'no indication annotated',
            'no target annotated'
        ])
    ]

    # Group by phase and calculate statistics
    forest = forest.groupby([f'{phase_type}num', phase_type]).agg(
        total_targets=pd.NamedAgg(
            column='target_status',
            aggfunc=lambda x: x.notna().sum()  # Count non-null target statuses
        ),
        number_targets_with_genetic_support=pd.NamedAgg(
            column='target_status',
            aggfunc=lambda x: (x == 'genetically supported target').sum()
        )
    ).reset_index()

    return forest


def deduplicate_phase_data(
    phase_data: pd.DataFrame,
    phase_type: str,
    score_column: Optional[str] = None,
    score_strategy: str = 'max'
) -> pd.DataFrame:
    """
    Deduplicate ti_uid pairs for a single phase.

    This implements complex deduplication logic to ensure each trait-indication
    pair (ti_uid) appears only once per phase, taking the best match.

    Deduplication steps:
    1. Group by ti_uid, take max(similarity) and max(phase_category)
    2. Filter to rows matching these max values
    3. Add 'highest_phase_with_known_outcome' column
    4. Sort by ti_uid and highest_phase_with_known_outcome (descending)
    5. Apply score_strategy to select final row per ti_uid

    Args:
        phase_data: DataFrame for a specific phase
        phase_type: One of 'ccat', 'hcat', 'acat'
        score_column: Optional column name containing association scores.
        score_strategy: Strategy for handling scores when deduplicating. One of:
            - 'max': Keep row with highest score (sort descending, keep first)
            - 'min': Keep row with lowest score (sort ascending, keep first)
            - 'avg': Average scores within each ti_uid group before deduplication

    Returns:
        Deduplicated DataFrame with one row per ti_uid

    Reference: relative_success.py:220-253

    Note:
        This logic is subtle and critical. The rationale is:
        - Take maximum similarity (best genetic match for this indication)
        - Take maximum phase (most advanced clinical trial stage)
        - Prioritize rows with known outcomes at advanced phases
        - This ensures we're using the best available data for each ti_uid
        If a score column is provided, the score_strategy determines how to
        handle multiple rows with the same ti_uid after filtering.

    Example:
        >>> data = pd.DataFrame({
        ...     'ti_uid': ['T1-I1', 'T1-I1', 'T1-I1'],
        ...     'similarity': [0.9, 0.8, 0.9],
        ...     'ccatnum': [2, 2, 1],
        ...     'ccat': ['Phase II', 'Phase II', 'Phase I'],
        ...     'succ_p_1': [1.0, 1.0, 1.0],
        ...     'succ_1_2': [1.0, 1.0, np.nan],
        ...     'succ_2_3': [1.0, np.nan, np.nan],
        ...     'succ_3_a': [np.nan, np.nan, np.nan]
        ... })
        >>> deduped = deduplicate_phase_data(data, 'ccat')
        >>> len(deduped)
        1
        >>> deduped['ccatnum'].iloc[0]
        2
        >>> deduped['similarity'].iloc[0]
        0.9
    """
    if score_strategy not in ('max', 'min', 'avg'):
        raise ValueError(f"score_strategy must be 'max', 'min', or 'avg', got '{score_strategy}'")
    
    # Step 1: Group by ti_uid and get max similarity and max phase category
    step1 = phase_data.groupby('ti_uid').agg(
        maxsim=pd.NamedAgg(column='similarity', aggfunc='max'),
        maxcat=pd.NamedAgg(column=f'{phase_type}num', aggfunc='max'),
    ).reset_index()

    # Step 2: Filter to rows matching max values
    # This keeps only the rows with the best similarity and most advanced phase
    phase_data = phase_data[
        (phase_data['ti_uid'].isin(step1['ti_uid'])) &
        (phase_data['similarity'].isin(step1['maxsim'].unique())) &
        (phase_data[f'{phase_type}num'].isin(step1['maxcat'].unique()))
    ]

    # Merge back to get all columns with the max values
    # Note: When merging with different left_on/right_on column names, pandas keeps
    # both sets of columns. We need to drop the right-side merge keys before renaming
    # to avoid creating duplicate column names.
    step2 = step1.merge(
        phase_data,
        left_on=['ti_uid', 'maxsim', 'maxcat'],
        right_on=['ti_uid', 'similarity', f'{phase_type}num'],
        how='left',
    )
    # Drop the right-side merge keys (they're duplicates of maxsim, maxcat)
    step2 = step2.drop(columns=['similarity', f'{phase_type}num'])
    # Now rename the left-side keys to their proper names
    step2 = step2.rename(columns={'maxsim': 'similarity', 'maxcat': f'{phase_type}num'})

    # Step 3: Prioritize rows with known outcome at most advanced phase
    # This helps break ties when multiple rows have same max similarity and phase
    # We want the row with the most complete outcome data
    step2['highest_phase_with_known_outcome'] = step2.apply(
        lambda row: 4 if pd.notna(row['succ_3_a']) else
                    3 if pd.notna(row['succ_2_3']) else
                    2 if pd.notna(row['succ_1_2']) else
                    1 if pd.notna(row['succ_p_1']) else None,
        axis=1
    )

    # Filter out rows with no known outcome
    # This fixes bug in original code where None values were being treated as True
    step2 = step2[step2['highest_phase_with_known_outcome'].notna()]

    # Step 4: Apply score strategy and deduplicate
    if score_column is not None:
        if score_strategy == 'avg':
            # For each ti_uid, find rows with the max highest_phase_with_known_outcome,
            # then average scores only among those rows
            max_phase_per_ti = step2.groupby('ti_uid')['highest_phase_with_known_outcome'].transform('max')
            step2_max_phase = step2[step2['highest_phase_with_known_outcome'] == max_phase_per_ti].copy()

            # Compute average score within each ti_uid (only among rows with max phase)
            avg_scores = step2_max_phase.groupby('ti_uid')[score_column].mean().reset_index()
            avg_scores = avg_scores.rename(columns={score_column: f'{score_column}_avg'})

            # Keep first row per ti_uid (they all have same max phase)
            step2 = step2_max_phase.drop_duplicates(subset=['ti_uid'], keep='first').reset_index(drop=True)

            # Merge averaged score back and replace original score column
            step2 = step2.merge(avg_scores, on='ti_uid', how='left')
            step2[score_column] = step2[f'{score_column}_avg']
            step2 = step2.drop(columns=[f'{score_column}_avg'])
        else:
            # 'max' or 'min' strategy: sort by score and keep first
            # For 'max': sort descending (ascending=False), keep first = highest
            # For 'min': sort ascending (ascending=True), keep first = lowest
            score_ascending = (score_strategy == 'min')
            step2 = step2.sort_values(
                by=['ti_uid', 'highest_phase_with_known_outcome', score_column],
                ascending=[True, False, score_ascending]
            )
            step2 = step2.drop_duplicates(subset=['ti_uid'], keep='first').reset_index(drop=True)
    else:
        # No score column - just sort by highest_phase_with_known_outcome
        step2 = step2.sort_values(
            by=['ti_uid', 'highest_phase_with_known_outcome'],
            ascending=[True, False]
        )
        step2 = step2.drop_duplicates(subset=['ti_uid'], keep='first').reset_index(drop=True)

    return step2

def generate_forest_data_for_phase_type(
    merged_df: pd.DataFrame,
    phase_type: str,
    merged_phases: Set[str]
) -> pd.DataFrame:
    """
    Generate complete forest data for one phase type (ccat/hcat/acat).

    Steps:
    1. For each unique phase in merged_phases:
       a. Filter to rows in that phase
       b. Deduplicate by ti_uid
       c. Generate forest data
    2. Concatenate all phase forest data
    3. Add any missing phases (0-4) with zero counts
    4. Sort by phase number

    Args:
        merged_df: Complete merged dataframe with all data
        phase_type: One of 'ccat', 'hcat', 'acat'
        merged_phases: Set of unique phase names found in data

    Returns:
        Forest data DataFrame for this phase type

    Reference: relative_success.py:213-270

    Note:
        Missing phases (0-4) are added with zero counts to ensure complete
        data for all phases. This is important for RS calculations which
        expect all 5 phases to be present.

    Example:
        >>> df = pd.DataFrame({
        ...     'ccat': ['Preclinical', 'Phase I'],
        ...     'ccatnum': [0, 1],
        ...     'ti_uid': ['T1-I1', 'T2-I2'],
        ...     'similarity': [0.9, 0.8],
        ...     'target_status': ['genetically supported target', 'unsupported target'],
        ...     'succ_p_1': [1.0, 1.0],
        ...     'succ_1_2': [np.nan, np.nan],
        ...     'succ_2_3': [np.nan, np.nan],
        ...     'succ_3_a': [np.nan, np.nan]
        ... })
        >>> forest = generate_forest_data_for_phase_type(df, 'ccat', {'Preclinical', 'Phase I'})
        >>> len(forest)  # Should have entries for all 5 phases (0-4)
        5
    """
    forest_data_list = []

    # Process each phase
    for phase in merged_phases:
        # Filter to this phase
        merged_phase = merged_df[merged_df[phase_type] == phase].copy()

        if merged_phase.empty:
            continue

        # Deduplicate by ti_uid
        deduplicated = deduplicate_phase_data(merged_phase, phase_type)

        # Generate forest statistics
        forest = get_forest_data(phase_type, deduplicated)

        forest_data_list.append(forest)

    # Concatenate all non-empty forest data
    try:
        forest_data = pd.concat(
            [df for df in forest_data_list if not df.empty and not df.isna().all(axis=None)],
            ignore_index=True
        )
    except ValueError:
        # Completely empty forest data
        print(f"Warning: Completely empty forest data for {phase_type}")
        # Create empty dataframe with correct structure
        # Use first non-empty dataframe's columns if available, otherwise create minimal structure
        if forest_data_list and not forest_data_list[0].empty:
            forest_data = forest_data_list[0].iloc[0:0]  # Empty df with same columns
        else:
            # Create minimal structure
            forest_data = pd.DataFrame(columns=[
                f'{phase_type}num',
                phase_type,
                'total_targets',
                'number_targets_with_genetic_support'
            ])

    # Add missing phases (0-4) with zero counts
    existing_phase_nums = [int(i) for i in forest_data[f'{phase_type}num'].unique()] if not forest_data.empty else []

    for phase_num in range(0, 5):
        if phase_num not in existing_phase_nums:
            print(f"Adding missing phase {phase_num} ({PHASE_MAPPING_REVERSE.get(phase_num, 'Unknown')}) to {phase_type} forest data")
            missing_row = pd.DataFrame({
                f'{phase_type}num': [phase_num],
                phase_type: [PHASE_MAPPING_REVERSE[phase_num]],
                'total_targets': [0],
                'number_targets_with_genetic_support': [0]
            })
            forest_data = pd.concat([forest_data, missing_row], ignore_index=True)

    # Sort by phase number
    forest_data = forest_data.sort_values(by=[f'{phase_type}num']).reset_index(drop=True)

    return forest_data
