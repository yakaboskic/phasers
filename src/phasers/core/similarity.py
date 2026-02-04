"""
Similarity calculation between indication and association MeSH IDs.

Handles:
- Single MeSH ID comparisons
- Multiple MeSH IDs (take maximum similarity)
- Missing/invalid MeSH IDs (return negative codes)

Reference: relative_success.py lines 45-80, other_clinical_success_metrics.py lines 37-78
"""

import pandas as pd
import numpy as np
import re
import ast
from typing import Union, List, Set
from tqdm import tqdm

# Initialize tqdm for pandas
tqdm.pandas()

# Module-level variables set by calling code
# These need to be set before calling get_similarity functions
similarity_matrix: pd.DataFrame = None
BAD_SIMILARITIES: Set[str] = set()


def get_similarity(
    indication_mesh_id: str,
    association_mesh_id: str
) -> float:
    """
    Get similarity score between two MeSH IDs.

    This function looks up the similarity score in the global similarity_matrix.
    It handles various edge cases where MeSH IDs may be missing or not found.

    Args:
        indication_mesh_id: MeSH ID from pharmaprojects indication
        association_mesh_id: MeSH ID from genetic association

    Returns:
        float: Similarity score or negative error code:
            - Positive (0-1): Actual similarity score from matrix
            - -1: Both IDs are NaN (no indication, no association)
            - -2: Indication ID is NaN (no pharmaprojects indication)
            - -3: Association ID is NaN (no genetic association found)
            - -4: MeSH ID not found in similarity matrix

    Reference: relative_success.py:65-80

    Note:
        These negative codes are used downstream for target_status assignment.
        The global similarity_matrix and BAD_SIMILARITIES must be set before calling.

    Example:
        >>> similarity_matrix = pd.DataFrame({'D003920': [1.0, 0.5]},
        ...                                   index=['D003920', 'D012345'])
        >>> get_similarity('D003920', 'D003920')
        1.0
        >>> get_similarity(np.nan, 'D003920')
        -2
    """
    global BAD_SIMILARITIES

    # Handle both MeSH IDs being NaN
    if pd.isna(indication_mesh_id) and pd.isna(association_mesh_id):
        # Pharma projects didn't have an indication for this target
        # AND no association was found for the target
        return -1

    # Handle indication MeSH ID being NaN
    elif pd.isna(indication_mesh_id):
        # Pharma projects didn't have an indication for this target
        return -2

    # Handle association MeSH ID being NaN
    elif pd.isna(association_mesh_id):
        # This indication didn't have an association in the data
        return -3

    # Look up similarity in matrix
    try:
        sim = similarity_matrix.loc[indication_mesh_id, association_mesh_id]
        return sim

    except KeyError as e:
        # MeSH ID not found in similarity matrix
        # Extract the problematic ID from error message and track it
        match = re.search(r"\(([^)]+)\)", str(e))
        if match:
            bad_id = match.group(1).strip("'")
            BAD_SIMILARITIES.add(bad_id)
        return -4


def get_similarity_wrapper(
    indication_mesh_id: str,
    association_mesh_id: Union[str, List[str]]
) -> float:
    """
    Wrapper to handle associations with multiple MeSH IDs.

    If association_mesh_id is a list, returns max similarity across all IDs.
    This handles cases where an association maps to multiple MeSH IDs - we
    want to take the best match.

    Handles both:
    - Python list objects: ['MESH:123', 'MESH:456']
    - String representations of lists: "['MESH:123', 'MESH:456']"

    Args:
        indication_mesh_id: MeSH ID from pharmaprojects indication
        association_mesh_id: MeSH ID or list of MeSH IDs from association

    Returns:
        float: Maximum similarity score if multiple IDs, or single similarity

    Reference: relative_success.py:45-63

    Example:
        >>> get_similarity_wrapper('D003920', ['D003920', 'D012345'])
        1.0  # Takes max similarity
        >>> get_similarity_wrapper('D003920', 'D003920')
        1.0
    """
    # Handle list of mesh IDs (Python list object)
    if isinstance(association_mesh_id, list):
        max_sim = -np.inf
        for mesh_id in association_mesh_id:
            sim = get_similarity(indication_mesh_id, mesh_id)
            if sim > max_sim:
                max_sim = sim
        return max_sim

    # Handle string representation of list
    elif isinstance(association_mesh_id, str) and association_mesh_id.startswith('['):
        # Parse string to list
        association_mesh_id = ast.literal_eval(association_mesh_id)
        max_sim = -np.inf
        for mesh_id in association_mesh_id:
            sim = get_similarity(indication_mesh_id, mesh_id)
            if sim > max_sim:
                max_sim = sim
        return max_sim

    # Handle single mesh ID
    else:
        return get_similarity(indication_mesh_id, association_mesh_id)


def calculate_all_similarities(
    merged_df: pd.DataFrame,
    similarity_matrix_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate similarity for all rows in merged dataframe.

    Adds 'similarity' column using progress_apply with tqdm to show progress.
    This can be a slow operation for large datasets, hence the progress bar.

    Args:
        merged_df: Merged pharmaprojects/associations dataframe
        similarity_matrix_df: Similarity matrix (wide format, index=Indication_MeSH_ID,
                              columns=Association_MeSH_ID)

    Returns:
        DataFrame with added 'similarity' column

    Reference: relative_success.py:514-522

    Side Effects:
        - Sets global similarity_matrix variable
        - Populates global BAD_SIMILARITIES set
        - Prints warnings if bad similarities found

    Note:
        Must call tqdm.pandas() before using this function (done at module level).

    Example:
        >>> merged = pd.DataFrame({'indication_mesh_id': ['D003920'],
        ...                        'mesh_id': ['D003920']})
        >>> sim_matrix = pd.DataFrame({'D003920': [1.0]}, index=['D003920'])
        >>> result = calculate_all_similarities(merged, sim_matrix)
        >>> result['similarity'].iloc[0]
        1.0
    """
    global similarity_matrix, BAD_SIMILARITIES

    # Reset BAD_SIMILARITIES for this calculation
    BAD_SIMILARITIES = set()

    # Set module-level similarity matrix
    similarity_matrix = similarity_matrix_df

    # Calculate similarities with progress bar
    print("Computing similarity between indication and association MeSH IDs")
    merged_df['similarity'] = merged_df.progress_apply(
        lambda x: get_similarity_wrapper(x['indication_mesh_id'], x['mesh_id']),
        axis=1
    )

    # Report any bad similarities found
    if len(BAD_SIMILARITIES) > 0:
        print(
            f"Found {len(BAD_SIMILARITIES)} MeSH IDs that were not in the "
            f"similarity matrix: {BAD_SIMILARITIES}"
        )

    return merged_df
