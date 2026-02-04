"""
NEW MODULE: Track clinical trial levels for trait-indication (ti_uid) pairs.

For each ti_uid, determine:
- Minimum clinical trial level reached
- Maximum clinical trial level reached

Rationale: Multiple drugs may target the same (trait, indication) pair.
We want to know the range of clinical success across all attempts.

This addresses the requirement from INITIAL.md task #1:
"A new file that is each ti_uid, which is a trait - indication (gene) pair
and what its lowest and highest obtained clinical trail level it got to."
"""

import pandas as pd
from .phase_mapping import PHASE_MAPPING_REVERSE


def extract_ti_uid_levels(
    merged_df: pd.DataFrame,
    phase_type: str = 'ccat'
) -> pd.DataFrame:
    """
    Extract min/max clinical trial levels for each ti_uid.

    Groups merged data by ti_uid and calculates:
    - min_phase_num: Minimum phase reached (numeric)
    - max_phase_num: Maximum phase reached (numeric)
    - min_phase_name: Minimum phase name
    - max_phase_name: Maximum phase name
    - n_attempts: Number of unique drugs targeting this ti_uid

    Args:
        merged_df: Merged pharmaprojects/associations dataframe
        phase_type: Which phase type to use ('ccat', 'hcat', or 'acat')

    Returns:
        DataFrame with columns:
        - ti_uid: Trait-indication unique identifier
        - min_phase_num, max_phase_num: Numeric phase codes
        - min_phase_name, max_phase_name: Phase names
        - n_attempts: Number of unique drug attempts

    Note:
        This is NEW functionality requested in INITIAL.md.
        Multiple drugs may target the same (trait, indication) pair.
        This shows the range of clinical success across all attempts.

    Example:
        >>> df = pd.DataFrame({
        ...     'ti_uid': ['T1-I1', 'T1-I1', 'T2-I2'],
        ...     'gene': ['G1', 'G1', 'G2'],
        ...     'ccatnum': [1, 3, 0]
        ... })
        >>> result = extract_ti_uid_levels(df, 'ccat')
        >>> result.loc[result['ti_uid'] == 'T1-I1', 'min_phase_num'].iloc[0]
        1
        >>> result.loc[result['ti_uid'] == 'T1-I1', 'max_phase_num'].iloc[0]
        3
        >>> result.loc[result['ti_uid'] == 'T1-I1', 'n_attempts'].iloc[0]
        2
    """
    phase_num_col = f'{phase_type}num'

    # Only keep rows that have a similarity of 0.8 or greater (these are the matches we want to track)
    merged_df = merged_df[merged_df['similarity'] >= 0.8]

    # Filter out rows with invalid phase data
    valid_df = merged_df[
        merged_df[phase_num_col].notna() &
        (merged_df[phase_num_col] >= 0)  # Exclude NaN phase (-1)
    ].copy()

    # Group by ti_uid and calculate min/max
    ti_tracking = valid_df.groupby('ti_uid').agg(
        min_phase_num=pd.NamedAgg(column=phase_num_col, aggfunc='min'),
        max_phase_num=pd.NamedAgg(column=phase_num_col, aggfunc='max'),
        n_attempts=pd.NamedAgg(column='gene', aggfunc='nunique')
    ).reset_index()
    
    # Map phase numbers to names
    ti_tracking['min_phase_name'] = ti_tracking['min_phase_num'].map(PHASE_MAPPING_REVERSE)
    ti_tracking['max_phase_name'] = ti_tracking['max_phase_num'].map(PHASE_MAPPING_REVERSE)

    # Convert phase numbers to int for cleaner output
    ti_tracking['min_phase_num'] = ti_tracking['min_phase_num'].astype(int)
    ti_tracking['max_phase_num'] = ti_tracking['max_phase_num'].astype(int)

    return ti_tracking


def generate_ti_tracking_file(
    merged_df: pd.DataFrame,
    output_file: str
) -> str:
    """
    Generate ti_uid tracking file for all phase types.

    Combines ti_uid tracking for ccat, hcat, acat with appropriate column prefixes.
    This shows the min/max clinical levels reached for each trait-indication pair
    across all three phase categorizations.

    Args:
        merged_df: Merged pharmaprojects/associations dataframe
        output_file: Path to output CSV file

    Returns:
        Path to generated file

    Note:
        The file includes separate columns for each phase type:
        - ccat_*: Combined category (all trials)
        - hcat_*: Historical category (past trials)
        - acat_*: Active category (ongoing trials)

    Example:
        >>> df = pd.DataFrame({
        ...     'ti_uid': ['T1-I1'],
        ...     'ccatnum': [2], 'hcatnum': [1], 'acatnum': [3]
        ... })
        >>> file_path = generate_ti_tracking_file(df, '/tmp/test.ti_tracking.csv')
        >>> os.path.exists(file_path)
        True
    """
    print("Generating ti_uid tracking file")

    # Extract tracking for each phase type
    ccat_tracking = extract_ti_uid_levels(merged_df, 'ccat')
    hcat_tracking = extract_ti_uid_levels(merged_df, 'hcat')
    acat_tracking = extract_ti_uid_levels(merged_df, 'acat')

    # Extract clinical areas from merged df
    areas = merged_df[['ti_uid', 'areas']].dropna().drop_duplicates()

    # Rename columns with phase type prefix
    ccat_tracking = ccat_tracking.rename(columns={
        'min_phase_num': 'ccat_min_phase_num',
        'max_phase_num': 'ccat_max_phase_num',
        'min_phase_name': 'ccat_min_phase_name',
        'max_phase_name': 'ccat_max_phase_name',
        'n_attempts': 'ccat_n_attempts'
    })

    hcat_tracking = hcat_tracking.rename(columns={
        'min_phase_num': 'hcat_min_phase_num',
        'max_phase_num': 'hcat_max_phase_num',
        'min_phase_name': 'hcat_min_phase_name',
        'max_phase_name': 'hcat_max_phase_name',
        'n_attempts': 'hcat_n_attempts'
    })

    acat_tracking = acat_tracking.rename(columns={
        'min_phase_num': 'acat_min_phase_num',
        'max_phase_num': 'acat_max_phase_num',
        'min_phase_name': 'acat_min_phase_name',
        'max_phase_name': 'acat_max_phase_name',
        'n_attempts': 'acat_n_attempts'
    })

    # Merge all phase types
    ti_tracking_combined = ccat_tracking.merge(
        hcat_tracking,
        on='ti_uid',
        how='outer'
    ).merge(
        acat_tracking,
        on='ti_uid',
        how='outer'
    )

    # Seperate out the ti_uid into trait and indication
    ti_tracking_combined['trait'] = ti_tracking_combined['ti_uid'].str.split('-').str[0]
    ti_tracking_combined['indication'] = ti_tracking_combined['ti_uid'].str.split('-').str[1]

    # Merge in clinical areas
    ti_tracking_combined = ti_tracking_combined.merge(
        areas,
        on='ti_uid',
        how='left'
    )

    # Write output file
    ti_tracking_combined.to_csv(output_file, index=False, sep='\t')

    print(f"Ti_uid tracking file written to {output_file}")
    print(f"  Total unique ti_uid pairs: {len(ti_tracking_combined)}")

    return output_file
