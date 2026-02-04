"""
Data loading and preprocessing for relative success calculations.

This module handles:
- Loading pharmaprojects, indications, associations, and similarity matrix files
- Applying association threshold filters
- Merging dataframes for RS analysis
- Adding NaN duplicate rows (minikel compatibility)

Reference: relative_success.py lines 410-577
"""

import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Dict
import logging


def parse_list_str(column_value):
    """
    Parse string representation of lists from CSV files.

    Handles association files where mesh_id column may be stored as string:
    "['MESH:123', 'MESH:456']"

    Args:
        column_value: Value from CSV, may be string representation of list

    Returns:
        List if parseable, otherwise original value

    Reference: relative_success.py:82-85, other_clinical_success_metrics.py:74-77

    Example:
        >>> parse_list_str("['D003920', 'D012345']")
        ['D003920', 'D012345']
        >>> parse_list_str('D003920')
        'D003920'
    """
    if isinstance(column_value, str) and column_value.startswith('[') and column_value.endswith(']'):
        return column_value[1:-1].replace("'", "").split(',')
    return column_value


def load_pharmaprojects(file_path: str) -> pd.DataFrame:
    """
    Load pharmaprojects clinical trial data.

    Expected columns:
    - gene: Target gene
    - indication_mesh_id: Indication MeSH ID
    - ccat, hcat, acat: Clinical trial phases (combined, historical, active)
    - succ_p_1, succ_1_2, succ_2_3, succ_3_a: Success indicators at each phase transition
    - ti_uid: Trait-indication unique identifier (optional, generated if not present)

    Args:
        file_path: Path to pharmaprojects TSV file

    Returns:
        DataFrame with pharmaprojects data

    Reference: relative_success.py:416-417

    Example:
        >>> pp = load_pharmaprojects('raw/pharmaprojects.tsv')
        >>> 'gene' in pp.columns
        True
    """
    print(f"Loading pharma projects data from {file_path}")
    pharmaprojects = pd.read_csv(file_path, sep='\t')
    return pharmaprojects


def load_indications(file_path: str) -> pd.DataFrame:
    """
    Load indications data mapping MeSH IDs to areas and genetic insight.

    Expected columns:
    - indication_mesh_id: MeSH ID
    - genetic_insight: Whether indication has genetic evidence ('none' or other values)
    - areas: Clinical area classifications (comma-separated string)

    Args:
        file_path: Path to indications TSV file

    Returns:
        DataFrame with indications data

    Reference: relative_success.py:419-420

    Example:
        >>> ind = load_indications('raw/indications.tsv')
        >>> 'genetic_insight' in ind.columns
        True
    """
    print(f"Loading indications data from {file_path}")
    indications = pd.read_csv(file_path, sep='\t')
    return indications


def load_associations(
    file_path: str,
    pigean_file_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Load genetic associations data.

    Expected columns vary by source, but normalized to include:
    - gene: Target gene
    - mesh_id: Association MeSH ID (may be list)
    - source: Association source (e.g., 'FALCON', 'MAGMA', 'pigean')
    - [source-specific score columns like pos_p, p_value, combined, etc.]

    Args:
        file_path: Path to main associations file
        pigean_file_path: Optional path to separate pigean associations file

    Returns:
        DataFrame with associations data (concatenated if pigean file provided)

    Reference: relative_success.py:422-433

    Note:
        If pigean_file_path is provided, both files are loaded and concatenated.
        The mesh_id column is parsed to handle string representations of lists.

    Example:
        >>> assoc = load_associations('raw/associations.tsv')
        >>> 'gene' in assoc.columns and 'mesh_id' in assoc.columns
        True
    """
    print(f"Loading associations data from {file_path}")

    # Load main associations file
    if file_path:
        associations = pd.read_csv(file_path, sep='\t')
    else:
        associations = pd.DataFrame()

    # Load and append pigean associations if provided
    if pigean_file_path:
        print(f"Loading pigean associations data from {pigean_file_path}")
        pigean_associations = pd.read_csv(
            pigean_file_path,
            sep='\t',
            converters={'mesh_id': parse_list_str}
        )
        associations = pd.concat([associations, pigean_associations], ignore_index=True)

    # Ensure we have associations data
    if file_path is None and pigean_file_path is None:
        raise ValueError("Either file_path or pigean_file_path must be provided")

    return associations


def load_similarity_matrix(
    file_path: str,
    format: str = 'long'
) -> pd.DataFrame:
    """
    Load and pivot similarity matrix.

    Input format options:
    - 'long': Three columns - Indication_MeSH_ID, Association_MeSH_ID, Combined_Similarity
    - 'matrix': Already pivoted with Indication_MeSH_ID as index

    Returns:
        Wide format DataFrame with Indication_MeSH_ID as index,
        Association_MeSH_ID as columns, Combined_Similarity as values

    Args:
        file_path: Path to similarity matrix file
        format: Format of input file ('long' or 'matrix')

    Returns:
        DataFrame in wide format for efficient lookup

    Reference: relative_success.py:436-446

    Note:
        Pivoting is required for efficient lookups during similarity calculations.
        The long format is more storage-efficient but needs to be pivoted.

    Example:
        >>> sim_matrix = load_similarity_matrix('raw/minikel.similarity.csv', 'long')
        >>> sim_matrix.loc['D003920', 'D003920']  # Self-similarity
        1.0
    """
    print(f"Loading similarity matrix from {file_path}")

    if format == 'long':
        # Load long format and pivot to wide
        similarity_matrix = pd.read_csv(file_path)
        similarity_matrix = similarity_matrix.pivot(
            index='Indication_MeSH_ID',
            columns='Association_MeSH_ID',
            values='Combined_Similarity'
        )
    elif format == 'matrix':
        # Already in wide format
        similarity_matrix = pd.read_csv(file_path, index_col='Indication_MeSH_ID')
    else:
        raise ValueError(f"Invalid similarity matrix format: {format}. Must be 'long' or 'matrix'")

    return similarity_matrix


def parse_association_threshold(
    threshold_str: str
) -> Tuple[str, str, str, float, Optional[float]]:
    """
    Parse association threshold filter specification.

    Format: "SOURCE:COLUMN:OPERATOR:VALUE1:VALUE2"
    Examples:
        - "FALCON:pos_p:>=:0.5:"
        - "pigean:combined:>=:3:"
        - "MAGMA:p_value:<=:1e-8:"

    Args:
        threshold_str: Threshold specification string

    Returns:
        Tuple of (source_name, score_column, operator, value1, value2)
        value2 is None if not provided

    Reference: relative_success.py:552-577

    Note:
        Supported operators: '<', '>', '<=', '>=', '<>' (outside range),
                            '><' (inside range), '==' (equal)

    Example:
        >>> parse_association_threshold("FALCON:pos_p:>=:0.5:")
        ('FALCON', 'pos_p', '>=', 0.5, None)
        >>> parse_association_threshold("gwas:log_p:<>:1:2")
        ('gwas', 'log_p', '<>', 1.0, 2.0)
    """
    parts = threshold_str.split(':')

    association_name = parts[0]
    association_score_column = parts[1]
    threshold_operation = parts[2]
    threshold_value_1 = float(parts[3])

    # Handle optional second value
    threshold_value_2 = None
    if len(parts) > 4 and parts[4]:
        threshold_value_2 = float(parts[4])

    return (
        association_name,
        association_score_column,
        threshold_operation,
        threshold_value_1,
        threshold_value_2
    )


def apply_association_filters(
    associations_data: Dict[str, pd.DataFrame],
    threshold_specs: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Apply association threshold filters to association dataframes.

    Filters are applied to source-specific dataframes before merging.
    Each threshold specification filters rows based on score column values.

    Args:
        associations_data: Dict mapping source names to DataFrames
        threshold_specs: List of threshold specification strings

    Returns:
        Dict of filtered DataFrames

    Reference: relative_success.py:552-577

    Note:
        Filtering happens BEFORE similarity calculation to reduce computation.
        Supports operators: '<', '>', '<=', '>=', '<>' (outside), '><' (inside), '=='

    Example:
        >>> data = {'FALCON': pd.DataFrame({'pos_p': [0.3, 0.7, 0.9]})}
        >>> filtered = apply_association_filters(data, ["FALCON:pos_p:>=:0.5:"])
        >>> len(filtered['FALCON'])
        2
    """
    print("Processing association thresholds")

    for threshold_str in threshold_specs:
        print(f"Processing association threshold: {threshold_str}")

        association_name, score_column, operation, value1, value2 = \
            parse_association_threshold(threshold_str)

        # Skip if this source isn't in our data
        if association_name not in associations_data:
            print(f"Warning: Association source '{association_name}' not found in data, skipping filter")
            continue

        df = associations_data[association_name]

        # Apply filter based on operation
        if operation == '<':
            associations_data[association_name] = df[df[score_column] < value1]
        elif operation == '>':
            associations_data[association_name] = df[df[score_column] > value1]
        elif operation == '<=':
            associations_data[association_name] = df[df[score_column] <= value1]
        elif operation == '>=':
            associations_data[association_name] = df[df[score_column] >= value1]
        elif operation == '<>':
            # Outside range: value < value1 OR value > value2
            associations_data[association_name] = df[
                (df[score_column] < value1) | (df[score_column] > value2)
            ]
        elif operation == '><':
            # Inside range: value1 < value < value2
            associations_data[association_name] = df[
                (df[score_column] > value1) & (df[score_column] < value2)
            ]
        elif operation == '==':
            associations_data[association_name] = df[df[score_column] == value1]
        else:
            raise ValueError(f"Invalid threshold operation: {operation}")

        print(f"  Filtered {association_name}: {len(df)} -> {len(associations_data[association_name])} rows")

    return associations_data


def apply_association_filters_to_merged_data(
    merged_data: pd.DataFrame,
    threshold_specs: List[str]
) -> pd.DataFrame:
    """
    Apply association threshold filters to merged data.

    Args:
        merged_data: Merged data DataFrame
        threshold_specs: List of threshold specification strings

    Returns:
        Merged data DataFrame with filters applied

    Reference: relative_success.py:552-577

    Note:
        Filtering happens AFTER similarity calculation to reduce computation.
        Supports operators: '<', '>', '<=', '>=', '<>' (outside), '><' (inside), '=='

    Example:
        >>> merged = apply_association_filters_to_merged_data(merged, ["FALCON:pos_p:>=:0.5:"])
    """
    print("Processing association thresholds")
    for threshold_str in threshold_specs:
        print(f"Processing association threshold: {threshold_str}")

        source_name, score_column, operation, value1, value2 = \
            parse_association_threshold(threshold_str)

        logging.info(f"Source name: {source_name}")
        logging.info(f'Sources in merged data: {merged_data["source"].unique()}')

        # Skip if this source isn't in our data
        if source_name not in merged_data['source'].unique():
            print(f"Warning: Association source '{source_name}' not found in data, skipping filter")
            continue

        # Split merged data into source vs all other sources
        source_data = merged_data[merged_data['source'] == source_name]
        og_len = len(source_data)
        other_data = merged_data[merged_data['source'] != source_name]

        # Apply filter based on operation
        if operation == '<':
            source_data = source_data[source_data[score_column] < value1]
        elif operation == '>':
            source_data = source_data[source_data[score_column] > value1]
        elif operation == '<=':
            source_data = source_data[source_data[score_column] <= value1]
        elif operation == '>=':
            source_data = source_data[source_data[score_column] >= value1]
        elif operation == '<>':
            # Outside range: value < value1 OR value > value2
            source_data = source_data[
                (source_data[score_column] < value1) | (source_data[score_column] > value2)
            ]
        elif operation == '><':
            # Inside range: value1 < value < value2
            source_data = source_data[
                (source_data[score_column] > value1) & (source_data[score_column] < value2)
            ]
        elif operation == '==':
            source_data = source_data[source_data[score_column] == value1]
        else:
            raise ValueError(f"Invalid threshold operation: {operation}")

        print(f"  Filtered {source_name}: {og_len} -> {len(source_data)} rows")
        merged_data = pd.concat([source_data, other_data], ignore_index=True)
    return merged_data


def merge_data(
    pharmaprojects: pd.DataFrame,
    associations: pd.DataFrame,
    indications: pd.DataFrame,
    no_genetic_insight: bool = False
) -> pd.DataFrame:
    """
    Merge pharmaprojects, associations, and indications data.

    Key steps:
    1. Merge pharmaprojects + associations on 'gene' (left join)
    2. Merge in genetic_insight and areas from indications
    3. Filter out empty genes/indication_mesh_ids
    4. Apply genetic insight filter if required
    5. Add duplicate pharmaprojects rows with NaN associations (minikel artifact)

    Args:
        pharmaprojects: Pharmaprojects clinical trial data
        associations: Genetic associations data
        indications: Indications mapping to areas and genetic insight
        no_genetic_insight: If False, filter to indications with genetic insight

    Returns:
        Merged DataFrame ready for similarity calculation

    Reference: relative_success.py:449-485

    Note:
        The NaN duplicate rows (step 5) are a minikel artifact but must be preserved
        for compatibility. They ensure all pharmaprojects targets are represented
        even without genetic associations.

    Example:
        >>> pp = pd.DataFrame({'gene': ['GENE1'], 'indication_mesh_id': ['D003920']})
        >>> assoc = pd.DataFrame({'gene': ['GENE1'], 'mesh_id': ['D003920']})
        >>> ind = pd.DataFrame({'indication_mesh_id': ['D003920'],
        ...                     'genetic_insight': 'gwas', 'areas': 'metabolic'})
        >>> merged = merge_data(pp, assoc, ind)
        >>> 'mesh_id' in merged.columns
        True
    """
    print("Merging pharma projects and associations data on gene (target)")

    # Step 1: Merge pharmaprojects + associations on gene
    merged = pd.merge(
        pharmaprojects,
        associations,
        on='gene',
        how='left',
        suffixes=('', '_association')
    )

    # Step 2: Merge in genetic_insight and areas from indications
    print("Merging in genetic insight column from indications data")
    merged = pd.merge(
        merged,
        indications[['indication_mesh_id', 'genetic_insight', 'areas']],
        on='indication_mesh_id',
        how='left'
    )

    # Step 3: Filter out rows with missing genes or indication_mesh_ids
    print("Filtering out rows with empty genes or indication_mesh_ids")
    merged = merged[
        (merged['gene'] != '') &
        (merged['indication_mesh_id'] != '') &
        (merged['gene'].notna()) &
        (merged['indication_mesh_id'].notna())
    ]

    # Step 4: Apply genetic insight filter if required
    if not no_genetic_insight:
        print('Applying genetic insight filter (keeping indications with genetic insight != "none")')
        merged = merged[
            merged['indication_mesh_id'].isin(
                indications[indications['genetic_insight'] != 'none']['indication_mesh_id']
            )
        ]
    else:
        print('Running without genetic insight filter (keeping all indications)')
        merged = merged[
            merged['indication_mesh_id'].isin(indications['indication_mesh_id'])
        ]

    # Step 5: Add duplicate pharmaprojects rows with NaN for associations columns
    # This is a minikel artifact - ensures all pharmaprojects targets are represented
    print("Adding duplicate of pharmaprojects rows with NaN for associations columns")
    pharma_nan = pharmaprojects.copy()

    # Set all association columns to NaN (except gene, which is the merge key)
    for col in associations.columns:
        if col != 'gene' and col != 'ti_uid':  # Keep gene and ti_uid intact
            pharma_nan[col] = np.nan

    # Merge in the indications data
    print("Merging indications into pharma_nan dataframe")
    pharma_nan = pd.merge(
        pharma_nan,
        indications[['indication_mesh_id', 'areas']],
        on='indication_mesh_id',
        how='left'
    )

    # Combine the two DataFrames
    print("Combining merged and pharma_nan DataFrames")
    merged = pd.concat([merged, pharma_nan], ignore_index=True)

    # Apply genetic insight filter again if required
    if not no_genetic_insight:
        merged = merged[
            merged['indication_mesh_id'].isin(
                indications[indications['genetic_insight'] != 'none']['indication_mesh_id']
            )
        ]

    return merged


def load_and_merge_all_data(
    pharmaprojects_file: str,
    indications_file: str,
    associations_file: str,
    similarity_matrix_file: str,
    pigean_file: Optional[str] = None,
    similarity_matrix_format: str = 'long',
    no_genetic_insight: bool = False,
    association_thresholds: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point: Load all data files and create merged dataframe.

    This is the primary function to call for loading and merging all RS data.
    It orchestrates loading, filtering, and merging of all input files.

    Args:
        pharmaprojects_file: Path to pharmaprojects TSV
        indications_file: Path to indications TSV
        associations_file: Path to associations TSV
        similarity_matrix_file: Path to similarity matrix CSV
        pigean_file: Optional path to separate pigean associations file
        similarity_matrix_format: Format of similarity matrix ('long' or 'matrix')
        no_genetic_insight: If True, don't filter by genetic insight
        association_thresholds: Optional list of threshold filter specifications

    Returns:
        Tuple of (merged_df, similarity_matrix_df)

    Example:
        >>> merged, sim_matrix = load_and_merge_all_data(
        ...     'raw/pharmaprojects.tsv',
        ...     'raw/indications.tsv',
        ...     'raw/associations.tsv',
        ...     'raw/minikel.similarity.csv'
        ... )
        >>> 'similarity' not in merged.columns  # Similarity not calculated yet
        True
    """
    # Load individual files
    pharmaprojects = load_pharmaprojects(pharmaprojects_file)
    indications = load_indications(indications_file)
    associations = load_associations(associations_file, pigean_file)
    similarity_matrix = load_similarity_matrix(similarity_matrix_file, similarity_matrix_format)

    # Apply association filtering if specified
    if association_thresholds:
        # First, separate associations by source
        association_sources = associations['source'].dropna().unique()
        associations_data = {}
        for source in association_sources:
            associations_data[source] = associations[associations['source'] == source].copy()

        # Apply filters
        associations_data = apply_association_filters(associations_data, association_thresholds)

        # Recombine filtered associations
        associations = pd.concat(list(associations_data.values()), ignore_index=True)

    # Merge all data
    merged = merge_data(pharmaprojects, associations, indications, no_genetic_insight)

    return merged, similarity_matrix
