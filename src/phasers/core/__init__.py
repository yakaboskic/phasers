"""
Phasers Core Module - Relative success calculation engine.

Provides the computational core for relative success calculations:
- data_loading: Load and merge input files
- similarity: Calculate MeSH ID similarities
- phase_mapping: Phase categorization utilities
- forest_generation: Generate forest plot data
- rs_calculation: Core RS calculations with confidence intervals
- continuous: Continuous RS curves with incremental updates
- ti_tracking: Track clinical levels for target-indication pairs
- output_formatting: Write result files
"""

# Data loading
from .data_loading import (
    load_pharmaprojects,
    load_indications,
    load_associations,
    load_similarity_matrix,
    parse_association_threshold,
    apply_association_filters,
    apply_association_filters_to_merged_data,
    merge_data,
    load_and_merge_all_data,
)

# Similarity calculation
from .similarity import (
    get_similarity,
    get_similarity_wrapper,
    calculate_all_similarities,
)

# Phase mapping
from .phase_mapping import (
    PHASE_MAPPING,
    PHASE_MAPPING_REVERSE,
    add_phase_numeric_columns,
    get_unique_phases,
    validate_phases,
)

# Forest generation
from .forest_generation import (
    get_forest_data,
    deduplicate_phase_data,
    generate_forest_data_for_phase_type,
)

# RS calculation
from .rs_calculation import (
    calculate_rs_for_forest,
    delta_method_rs_ci,
    calculate_overall_rs,
)

# Continuous RS (incremental ROC-style)
from .continuous import (
    extract_deduplicated_records,
    build_rs_curve_incremental,
    calculate_summary_metrics,
    calculate_continuous_rs_incremental,
    write_continuous_rs_points,
)

# Ti_uid tracking
from .ti_tracking import (
    extract_ti_uid_levels,
    generate_ti_tracking_file,
)

# Output formatting
from .output_formatting import (
    write_forest_files,
    write_rs_summary,
    write_rs_summary_multi_source,
    write_contributions_file,
    write_merged_data,
    write_consolidated_forest_json,
)

__all__ = [
    # Data loading
    'load_pharmaprojects',
    'load_indications',
    'load_associations',
    'load_similarity_matrix',
    'parse_association_threshold',
    'apply_association_filters',
    'apply_association_filters_to_merged_data',
    'merge_data',
    'load_and_merge_all_data',
    # Similarity
    'get_similarity',
    'get_similarity_wrapper',
    'calculate_all_similarities',
    # Phase mapping
    'PHASE_MAPPING',
    'PHASE_MAPPING_REVERSE',
    'add_phase_numeric_columns',
    'get_unique_phases',
    'validate_phases',
    # Forest generation
    'get_forest_data',
    'deduplicate_phase_data',
    'generate_forest_data_for_phase_type',
    # RS calculation
    'calculate_rs_for_forest',
    'delta_method_rs_ci',
    'calculate_overall_rs',
    # Continuous RS
    'extract_deduplicated_records',
    'build_rs_curve_incremental',
    'calculate_summary_metrics',
    'calculate_continuous_rs_incremental',
    'write_continuous_rs_points',
    # Ti tracking
    'extract_ti_uid_levels',
    'generate_ti_tracking_file',
    # Output formatting
    'write_forest_files',
    'write_rs_summary',
    'write_rs_summary_multi_source',
    'write_contributions_file',
    'write_merged_data',
    'write_consolidated_forest_json',
]
