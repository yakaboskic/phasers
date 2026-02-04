"""
Phasers - Clinical Phase Relative Success Calculator

Calculate Relative Success (RS) metrics to evaluate how well genetic evidence
predicts drug target success across clinical trial phases.

Quick Start:
    >>> from phasers import calculate_continuous_rs, load_and_merge_data
    >>>
    >>> # Load and merge data
    >>> merged_df, sim_matrix = load_and_merge_data(
    ...     pharmaprojects_file='pharmaprojects.tsv',
    ...     indications_file='indications.tsv',
    ...     associations_file='associations.tsv',
    ...     similarity_matrix_file='similarity.csv'
    ... )
    >>>
    >>> # Calculate RS curve
    >>> results = calculate_continuous_rs(
    ...     merged_df=merged_df,
    ...     score_column='combined',
    ...     source_name='pigean'
    ... )

CLI Usage:
    $ phasers continuous --help
    $ phasers plot --help
"""

__version__ = "0.1.0"

# Import core functionality for convenient access
from .core import (
    # Data loading
    load_and_merge_all_data as load_and_merge_data,
    load_pharmaprojects,
    load_indications,
    load_associations,
    load_similarity_matrix,
    apply_association_filters,
    apply_association_filters_to_merged_data,
    # Similarity
    calculate_all_similarities,
    # Phase mapping
    PHASE_MAPPING,
    PHASE_MAPPING_REVERSE,
    add_phase_numeric_columns,
    # Forest generation
    generate_forest_data_for_phase_type,
    deduplicate_phase_data,
    # RS calculation
    calculate_rs_for_forest,
    delta_method_rs_ci,
    calculate_overall_rs,
    # Continuous RS
    calculate_continuous_rs_incremental as calculate_continuous_rs,
    build_rs_curve_incremental,
    calculate_summary_metrics,
    write_continuous_rs_points,
    # Ti tracking
    generate_ti_tracking_file,
    # Output
    write_forest_files,
    write_rs_summary,
    write_merged_data,
)

# Plotting
from .plotting import save_figures

__all__ = [
    # Version
    '__version__',
    # Data loading
    'load_and_merge_data',
    'load_pharmaprojects',
    'load_indications',
    'load_associations',
    'load_similarity_matrix',
    'apply_association_filters',
    'apply_association_filters_to_merged_data',
    # Similarity
    'calculate_all_similarities',
    # Phase mapping
    'PHASE_MAPPING',
    'PHASE_MAPPING_REVERSE',
    'add_phase_numeric_columns',
    # Forest generation
    'generate_forest_data_for_phase_type',
    'deduplicate_phase_data',
    # RS calculation
    'calculate_rs_for_forest',
    'delta_method_rs_ci',
    'calculate_overall_rs',
    # Continuous RS
    'calculate_continuous_rs',
    'build_rs_curve_incremental',
    'calculate_summary_metrics',
    'write_continuous_rs_points',
    # Ti tracking
    'generate_ti_tracking_file',
    # Output
    'write_forest_files',
    'write_rs_summary',
    'write_merged_data',
    # Plotting
    'save_figures',
]
