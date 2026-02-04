"""
Phase mapping and categorization utilities.

Handles conversion between phase names and numeric codes,
and extraction of unique phases from data.

Reference: relative_success.py lines 41-42, 500-511
"""

import pandas as pd
import numpy as np
from typing import Dict, Set


# Constants from relative_success.py:41-42
# Maps phase names to numeric codes for easier manipulation
PHASE_MAPPING: Dict = {
    np.nan: -1,
    'Preclinical': 0,
    'Phase I': 1,
    'Phase II': 2,
    'Phase III': 3,
    'Launched': 4
}

# Reverse mapping for converting numeric codes back to names
PHASE_MAPPING_REVERSE: Dict[int, str] = {
    0: 'Preclinical',
    1: 'Phase I',
    2: 'Phase II',
    3: 'Phase III',
    4: 'Launched'
}


def add_phase_numeric_columns(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add numeric phase columns (ccatnum, hcatnum, acatnum) to merged dataframe.

    Maps phase names to numeric codes for easier manipulation and calculations.
    The numeric codes allow for easy comparison and ordering of phases.

    Args:
        merged_df: DataFrame with ccat, hcat, acat columns containing phase names

    Returns:
        DataFrame with added ccatnum, hcatnum, acatnum columns

    Reference: relative_success.py:509-511

    Example:
        >>> df = pd.DataFrame({'ccat': ['Preclinical', 'Phase I', 'Launched']})
        >>> df = add_phase_numeric_columns(df)
        >>> df['ccatnum'].tolist()
        [0, 1, 4]
    """
    merged_df['ccatnum'] = merged_df['ccat'].map(PHASE_MAPPING)
    merged_df['hcatnum'] = merged_df['hcat'].map(PHASE_MAPPING)
    merged_df['acatnum'] = merged_df['acat'].map(PHASE_MAPPING)
    return merged_df


def get_unique_phases(merged_df: pd.DataFrame) -> Set[str]:
    """
    Get all unique phase values across ccat, hcat, acat columns.

    Used to validate all phases are covered in PHASE_MAPPING and to
    understand what phases exist in the data.

    Args:
        merged_df: DataFrame with ccat, hcat, acat columns

    Returns:
        Set of unique phase names found across all phase columns

    Reference: relative_success.py:500-508

    Example:
        >>> df = pd.DataFrame({'ccat': ['Preclinical', 'Phase I'],
        ...                    'hcat': ['Phase II'],
        ...                    'acat': ['Launched']})
        >>> phases = get_unique_phases(df)
        >>> 'Launched' in phases
        True
    """
    combined_phases = merged_df['ccat'].dropna().unique()
    historical_phases = merged_df['hcat'].dropna().unique()
    active_phases = merged_df['acat'].dropna().unique()

    # Union of all phase sets
    merged_phases = set.union(
        set(combined_phases),
        set(historical_phases),
        set(active_phases)
    )

    return merged_phases


def validate_phases(phases: Set[str]) -> None:
    """
    Validate that all phases in data exist in PHASE_MAPPING.

    Raises ValueError if an unknown phase is found. This prevents
    runtime errors when attempting to map phases to numeric codes.

    Args:
        phases: Set of phase names to validate

    Raises:
        ValueError: If any phase is not found in PHASE_MAPPING

    Reference: relative_success.py:505-508

    Example:
        >>> validate_phases({'Preclinical', 'Phase I'})  # No error
        >>> validate_phases({'Preclinical', 'Invalid Phase'})  # Raises ValueError
        Traceback (most recent call last):
        ...
        ValueError: Phase Invalid Phase not found in the hard-coded phase mapping
    """
    for phase in phases:
        if phase not in PHASE_MAPPING:
            raise ValueError(
                f"Phase {phase} not found in the hard-coded phase mapping. "
                f"Valid phases are: {list(PHASE_MAPPING.keys())}"
            )
