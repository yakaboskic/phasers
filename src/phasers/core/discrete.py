"""
Discrete relative success calculation logic.

Moved from run_workflow.py to phasers core so it can be reused
by both run-level and instance-level workflows.

Provides:
- annotate_target_status: Annotate target_status based on similarity threshold
- generate_baseline_associations: Create trait baseline associations
- calculate_discrete_rs: Calculate RS for all source/subsource/area combos
"""

import ast
import logging
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from tqdm import tqdm

from .forest_generation import generate_forest_data_for_phase_type
from .phase_mapping import add_phase_numeric_columns, get_unique_phases, validate_phases
from .rs_calculation import calculate_overall_rs, calculate_rs_for_forest


def annotate_target_status(
    merged_df: pd.DataFrame,
    similarity_threshold: float = 0.8,
) -> pd.DataFrame:
    """
    Annotate target_status based on similarity threshold.

    Target status categories:
    - 'genetically supported target': similarity >= threshold
    - 'unsupported target': similarity < threshold
    - 'indication lacks genetic insight': genetic_insight == 'none' (disabled)
    - 'no target annotated': gene is missing/empty
    - 'no indication annotated': indication_mesh_id is missing/empty

    Args:
        merged_df: Merged dataframe with similarity column.
        similarity_threshold: Threshold for genetic support.

    Returns:
        DataFrame with target_status column added.
    """
    logging.info(
        f"Annotating target status (similarity threshold: {similarity_threshold})"
    )

    merged_df['target_status'] = None

    merged_df.loc[
        merged_df['similarity'] >= similarity_threshold,
        'target_status',
    ] = 'genetically supported target'

    merged_df.loc[
        merged_df['similarity'] < similarity_threshold,
        'target_status',
    ] = 'unsupported target'

    # genetic_insight filtering is typically disabled
    merged_df.loc[
        False & (merged_df['genetic_insight'] == 'none'),
        'target_status',
    ] = 'indication lacks genetic insight'

    merged_df.loc[
        merged_df['gene'].isna() | (merged_df['gene'] == ''),
        'target_status',
    ] = 'no target annotated'

    merged_df.loc[
        merged_df['indication_mesh_id'].isna()
        | (merged_df['indication_mesh_id'] == ''),
        'target_status',
    ] = 'no indication annotated'

    status_counts = merged_df['target_status'].value_counts()
    logging.info("Target status distribution:")
    for status, count in status_counts.items():
        logging.info(f"  {status}: {count}")

    return merged_df


def generate_baseline_associations(
    merged_df: pd.DataFrame,
    pharmaprojects_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Generate trait baseline associations from merged data and pharmaprojects.

    For each (mesh_id, source, subsource) in the merged data, creates
    associations for every gene in pharmaprojects that targets that mesh_id.

    Args:
        merged_df: Merged dataframe (needs mesh_id, source, subsource columns).
        pharmaprojects_df: Pharmaprojects dataframe (needs indication_mesh_id, gene).

    Returns:
        DataFrame with columns: gene, mesh_id, source, subsource.
    """
    logging.info("Generating trait baseline associations")

    baseline_associations = merged_df[
        ['mesh_id', 'source', 'subsource']
    ].drop_duplicates()

    trait_baseline_associations = []
    for _, row in baseline_associations.iterrows():
        mesh_id = row['mesh_id']
        if pd.isna(mesh_id):
            continue
        try:
            mesh_id = ast.literal_eval(mesh_id)
        except Exception:
            mesh_id = str(mesh_id)
        source = row['source']
        subsource = row['subsource']
        if isinstance(mesh_id, list):
            for mesh_id_item in mesh_id:
                genes = pharmaprojects_df[
                    pharmaprojects_df['indication_mesh_id'] == mesh_id_item
                ]['gene'].unique()
                _df = pd.DataFrame(
                    {
                        'gene': genes,
                        'mesh_id': mesh_id_item,
                        'source': source,
                        'subsource': subsource,
                    }
                )
                trait_baseline_associations.append(_df)
        else:
            genes = pharmaprojects_df[
                pharmaprojects_df['indication_mesh_id'] == mesh_id
            ]['gene'].unique()
            _df = pd.DataFrame(
                {
                    'gene': genes,
                    'mesh_id': mesh_id,
                    'source': source,
                    'subsource': subsource,
                }
            )
            trait_baseline_associations.append(_df)

    if not trait_baseline_associations:
        return pd.DataFrame(columns=['gene', 'mesh_id', 'source', 'subsource'])

    result = pd.concat(trait_baseline_associations)
    result.drop_duplicates(inplace=True)
    logging.info(f"Generated {len(result)} baseline associations")
    return result


def calculate_discrete_rs(
    merged_df: pd.DataFrame,
    similarity_threshold: float = 0.8,
    use_subsources: bool = False,
    within_areas: bool = False,
) -> Tuple[
    Dict[str, Dict[Tuple[str, str], float]],
    Dict[str, Dict[str, pd.DataFrame]],
]:
    """
    Calculate discrete relative success for all source/subsource/area combos.

    Source naming convention: {source|all}__{subsource|all}__{area|all}

    Args:
        merged_df: Merged dataframe with similarity and target_status columns.
        similarity_threshold: Threshold for genetic support annotation.
        use_subsources: If True, calculate RS per subsource.
        within_areas: If True, calculate RS per clinical area.

    Returns:
        Tuple of (rs_results, forest_results):
        - rs_results: Dict[source_key, Dict[(phase_type, metric), value]]
        - forest_results: Dict[source_key, Dict[phase_type, DataFrame]]
    """
    logging.info("Calculating discrete relative success for all sources")

    # Annotate target status
    merged_df = annotate_target_status(merged_df, similarity_threshold)

    # Get unique sources (exclude intOGen — somatic, not germline)
    association_sources = merged_df['source'].dropna().unique()
    association_sources = [s for s in association_sources if s != 'intOGen']
    logging.info(f"Found association sources: {association_sources}")

    # Collect subsources if requested
    subsources: Optional[Dict[str, List[str]]] = None
    if use_subsources:
        _subsources = merged_df[['source', 'subsource']].dropna().drop_duplicates()
        subsources = {}
        for source, subsource in _subsources.values:
            if source not in subsources:
                subsources[source] = []
            subsources[source].append(subsource)
        logging.info(f"Found subsources: {subsources}")

    # Collect unique areas
    merged_df['areas'] = merged_df['areas'].fillna('None')
    areas: Set[str] = {
        y.strip() for x in merged_df['areas'].unique() for y in x.split(',')
    }
    if within_areas:
        logging.info(f"Found {len(areas)} clinical areas")

    # Phase validation
    merged_phases = get_unique_phases(merged_df)
    validate_phases(merged_phases)
    merged_df = add_phase_numeric_columns(merged_df)

    rs_results: Dict[str, Dict[Tuple[str, str], float]] = {}
    forest_results: Dict[str, Dict[str, pd.DataFrame]] = {}

    def _calc(task_df: pd.DataFrame, task_name: str) -> None:
        """Calculate RS for a data subset and store results."""
        logging.debug(f"  Calculating RS for: {task_name}")
        forest_dict: Dict[str, pd.DataFrame] = {}
        for phase_type in ['ccat', 'hcat', 'acat']:
            forest = generate_forest_data_for_phase_type(
                task_df, phase_type, merged_phases
            )
            forest = calculate_rs_for_forest(forest, phase_type)
            forest_dict[phase_type] = forest
        rs = calculate_overall_rs(forest_dict)
        rs_results[task_name] = rs
        forest_results[task_name] = forest_dict

    # All data
    _calc(merged_df, 'all__all__all')

    # Per source
    for source in tqdm(association_sources, desc='Calculating RS per source'):
        source_df = merged_df[
            (merged_df['source'] == source) | (merged_df['source'].isna())
        ]
        source_lower = source.lower()

        if use_subsources and subsources and source in subsources:
            for subsource in subsources[source]:
                subsource_df = source_df[
                    (source_df['subsource'] == subsource)
                    | (source_df['subsource'].isna())
                ]
                subsource_lower = subsource.lower()
                _calc(subsource_df, f'{source_lower}__{subsource_lower}__all')

            _calc(source_df, f'{source_lower}__all__all')
        else:
            _calc(source_df, f'{source_lower}__all__all')

    # Per area
    if within_areas:
        for area in tqdm(areas, desc='Calculating RS within areas'):
            if area == 'None':
                continue
            area_df = merged_df[merged_df['areas'].str.contains(area)]
            area_lower = area.lower()

            _calc(area_df, f'all__all__{area_lower}')

            for source in association_sources:
                source_area_df = area_df[
                    (area_df['source'] == source) | (area_df['source'].isna())
                ]
                source_lower = source.lower()

                if use_subsources and subsources and source in subsources:
                    for subsource in subsources[source]:
                        subsource_area_df = source_area_df[
                            (source_area_df['subsource'] == subsource)
                            | (source_area_df['subsource'].isna())
                        ]
                        subsource_lower = subsource.lower()
                        _calc(
                            subsource_area_df,
                            f'{source_lower}__{subsource_lower}__{area_lower}',
                        )
                    _calc(source_area_df, f'{source_lower}__all__{area_lower}')
                else:
                    _calc(source_area_df, f'{source_lower}__all__{area_lower}')

    return rs_results, forest_results
