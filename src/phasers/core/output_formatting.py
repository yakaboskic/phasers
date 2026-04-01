"""
Output formatting and file writing for RS results.

Handles:
- Writing forest CSV files
- Writing main RS summary CSV
- Formatting headers (flatten or multiindex)
- Creating log files

Reference: relative_success.py lines 386, 643-648
"""

import csv
import json
import logging
import os
from typing import Dict, List, Tuple

import pandas as pd


def write_forest_files(
    forest_dict: Dict[str, pd.DataFrame],
    output_dir: str,
    run_name: str,
    task_name: str = ''
) -> None:
    """
    Write forest data CSV files for each phase type.

    Creates files:
    - {output_dir}/{run_name}.relative_success.{task_name}.ccat.csv
    - {output_dir}/{run_name}.relative_success.{task_name}.hcat.csv
    - {output_dir}/{run_name}.relative_success.{task_name}.acat.csv

    Args:
        forest_dict: Dict mapping phase_type ('ccat', 'hcat', 'acat') to DataFrames
        output_dir: Directory to write files
        run_name: Prefix for output filenames
        task_name: Optional task identifier to include in filename

    Returns:
        None

    Reference: relative_success.py:386

    Example:
        >>> forests = {
        ...     'ccat': pd.DataFrame({'ccatnum': [0, 1], 'ccat': ['Preclinical', 'Phase I']})
        ... }
        >>> write_forest_files(forests, '/tmp', 'test', 'all')
        # Creates /tmp/test.relative_success.all.ccat.csv
    """
    for phase_type, forest_df in forest_dict.items():
        # Build filename
        if task_name:
            filename = f'{run_name}.relative_success.{task_name}.{phase_type}.csv'
        else:
            filename = f'{run_name}.relative_success.all.{phase_type}.csv'

        filepath = os.path.join(output_dir, filename)

        # Write CSV
        forest_df.to_csv(filepath, index=False)
        logging.info(f"  Wrote forest file: {filepath}")


def write_rs_summary(
    rs_dict: Dict[Tuple[str, str], float],
    output_dir: str,
    run_name: str,
    flatten_headers: bool = False
) -> None:
    """
    Write main relative success summary CSV.

    Converts rs_dict to DataFrame with MultiIndex columns or flattened headers.

    Format:
    - Rows: source names (e.g., 'all', 'FALCON', 'pigean')
    - Columns: (phase_type, metric_type) e.g., ('ccat', 'exact'), ('ccat', 'lower_ci')

    If flatten_headers=True: Join tuples to 'ccat-exact', 'ccat-lower_ci'

    Args:
        rs_dict: Dict with keys (phase_type, metric_type) and float values
        output_dir: Directory to write file
        run_name: Prefix for output filename
        flatten_headers: If True, flatten MultiIndex columns

    Returns:
        None

    Reference: relative_success.py:643-648

    Note:
        The rs_dict should have structure like:
        {
            ('ccat', 'exact'): 1.5,
            ('ccat', 'lower_ci'): 1.2,
            ('ccat', 'upper_ci'): 1.8,
            ...
        }

    Example:
        >>> rs_dict = {('ccat', 'exact'): 1.5, ('ccat', 'lower_ci'): 1.2}
        >>> write_rs_summary(rs_dict, '/tmp', 'test', flatten_headers=True)
        # Creates /tmp/test.relative_success.csv with 'ccat-exact' column
    """
    # Convert dict to DataFrame
    # Each key is (phase_type, metric_type), value is the RS value
    # We want to create a DataFrame where columns are MultiIndex

    # Group by row (we'll call it 'all' for now, more sophisticated logic in workflow)
    # For now, just create a single-row dataframe
    rs_df = pd.DataFrame([rs_dict]).T.reset_index()
    rs_df.columns = ['index', 'relative_success']

    # Split the index tuple into two columns
    rs_df[['phase_type', 'metric_type']] = pd.DataFrame(
        rs_df['index'].tolist(),
        index=rs_df.index
    )
    rs_df = rs_df.drop(columns=['index'])

    # Pivot to wide format
    rs_df_wide = rs_df.pivot(
        columns=['phase_type', 'metric_type'],
        values='relative_success'
    )

    # Add a row index (source name)
    # This will be set by the calling code, for now use 'all'
    rs_df_wide.index = ['all']
    rs_df_wide.index.name = 'source'

    # Flatten headers if requested
    if flatten_headers:
        # Join MultiIndex column tuples with '-'
        rs_df_wide.columns = ['-'.join(col).strip() if col[1] else col[0]
                               for col in rs_df_wide.columns.values]

    # Write CSV
    filepath = os.path.join(output_dir, f'{run_name}.relative_success.csv')
    rs_df_wide.to_csv(filepath, index=True)
    logging.info(f'Saved relative success to {filepath}')


def write_rs_summary_multi_source(
    rs_results: Dict[str, Dict[Tuple[str, str], float]],
    output_file: str,
    flatten_headers: bool = False
) -> None:
    """
    Write RS summary CSV with multiple sources/tasks.

    This is the more complete version that handles multiple sources
    (e.g., 'all', 'FALCON', 'pigean', etc.)

    Args:
        rs_results: Dict mapping source name to rs_dict
                    e.g., {'all': {('ccat', 'exact'): 1.5, ...},
                           'FALCON': {('ccat', 'exact'): 1.3, ...}}
        output_file: Path to output CSV file
        flatten_headers: If True, flatten MultiIndex columns

    Returns:
        None

    Example:
        >>> rs_results = {
        ...     'all': {('ccat', 'exact'): 1.5},
        ...     'FALCON': {('ccat', 'exact'): 1.3}
        ... }
        >>> write_rs_summary_multi_source(rs_results, '/tmp/test.relative_success.csv')
        # Creates /tmp/test.relative_success.csv with rows for 'all' and 'FALCON'
    """
    # Convert nested dict to DataFrame
    rs_df = pd.DataFrame(rs_results).T
    rs_df.index.name = 'source'

    # Flatten headers if requested
    if flatten_headers:
        # Join MultiIndex column tuples with '-'
        rs_df.columns = ['-'.join(str(c) for c in col).strip() if isinstance(col, tuple) else col
                         for col in rs_df.columns.values]

    # Write CSV
    rs_df.to_csv(output_file, index=True)
    logging.info(f'Saved relative success to {output_file}')


def write_contributions_file(
    contributions: Dict,
    output_dir: str,
    run_name: str
) -> None:
    """
    Write contributions CSV (if not disabled).

    Contributions show how each target-indication pair contributes to
    the overall relative success calculation.

    Args:
        contributions: Dict of contribution DataFrames
        output_dir: Directory to write file
        run_name: Prefix for output filename

    Returns:
        None

    Reference: relative_success.py:634-641

    Note:
        This functionality can be skipped if --no-contributions flag is set.
        Contributions are computationally expensive and not always needed.

    Example:
        >>> contribs = {
        ...     'all': {
        ...         'ccat': pd.DataFrame({'ti_uid': ['T1-I1'], 'contribution': [0.1]})
        ...     }
        ... }
        >>> write_contributions_file(contribs, '/tmp', 'test')
        # Creates /tmp/test.contributions.csv
    """
    # Concatenate all contribution dataframes
    all_contributions = []
    for source, contrib_dict in contributions.items():
        for phase_type, contrib_df in contrib_dict.items():
            contrib_df = contrib_df.copy()
            contrib_df['source'] = source
            contrib_df['phase_type'] = phase_type
            all_contributions.append(contrib_df)

    if all_contributions:
        contributions_df = pd.concat(all_contributions, ignore_index=True)

        # Write CSV
        filepath = os.path.join(output_dir, f'{run_name}.contributions.csv')
        contributions_df.to_csv(filepath, index=False)
        logging.info(f'Saved contributions to {filepath}')
    else:
        logging.info('No contributions to write')


def write_merged_data(
    merged_df: pd.DataFrame,
    output_file: str,
    compress: bool = True
) -> None:
    """
    Write merged data file (output of merge stage).

    Args:
        merged_df: Merged pharmaprojects/associations/similarity data
        output_file: Output filepath
        compress: If True, use gzip compression

    Returns:
        None

    Note:
        This is typically a large file, so compression is recommended.
        The file is tab-separated for compatibility with original scripts.

    Example:
        >>> merged = pd.DataFrame({'gene': ['G1'], 'similarity': [0.9]})
        >>> write_merged_data(merged, '/tmp/merged.tsv.gz', compress=True)
        # Creates /tmp/merged.tsv.gz
    """
    logging.info(f"Writing merged data to {output_file}")

    compression = 'gzip' if compress else None
    merged_df.to_csv(output_file, sep='\t', index=False, compression=compression)

    logging.info(f"Merged data saved to {output_file}")
    logging.info(f"  Rows: {len(merged_df)}")
    logging.info(f"  Columns: {len(merged_df.columns)}")


def write_consolidated_forest_json(
    all_forest_results: Dict[str, Dict[str, pd.DataFrame]],
    ccat_file: str,
    hcat_file: str,
    acat_file: str
) -> None:
    """
    Write consolidated forest data as JSON files (one per phase type).

    Instead of writing separate CSV files for each source/subsource/area,
    this writes one JSON file per phase type containing all forest data.

    Args:
        all_forest_results: Dict mapping source_name to forest_dict
                           e.g., {'all': {'ccat': df, 'hcat': df, 'acat': df},
                                  'otg': {'ccat': df, ...}, ...}
        ccat_file: Path to output ccat JSON file
        hcat_file: Path to output hcat JSON file
        acat_file: Path to output acat JSON file

    JSON structure:
        {
          "all": [{"ccatnum": 0, "ccat": "Preclinical", ...}, ...],
          "otg": [...],
          "otg_finngen": [...]
        }

    Example:
        >>> forests = {
        ...     'all': {'ccat': pd.DataFrame({'ccatnum': [0], 'ccat': ['Preclinical']})},
        ...     'otg': {'ccat': pd.DataFrame({'ccatnum': [0], 'ccat': ['Preclinical']})}
        ... }
        >>> write_consolidated_forest_json(forests, '/tmp/test.forest.ccat.json',
        ...                                '/tmp/test.forest.hcat.json',
        ...                                '/tmp/test.forest.acat.json')
    """
    logging.info("Writing consolidated forest JSON files")

    # Map phase types to output files
    phase_files = {
        'ccat': ccat_file,
        'hcat': hcat_file,
        'acat': acat_file
    }

    for phase_type, output_file in phase_files.items():
        # Collect all forest data for this phase type across all sources
        phase_data = {}

        for source_name, forest_dict in all_forest_results.items():
            if phase_type in forest_dict:
                # Convert DataFrame to records format (list of dicts)
                phase_data[source_name] = forest_dict[phase_type].to_dict('records')

        # Write JSON file
        with open(output_file, 'w') as f:
            json.dump(phase_data, f, indent=2)

        logging.info(f"  Wrote {output_file} ({len(phase_data)} sources)")

    logging.info("Consolidated forest JSON files complete")


def _parse_source_key(key: str) -> Tuple[str, str, str]:
    """
    Split a source__subsource__area key into its components.

    Args:
        key: String in format "source__subsource__area".

    Returns:
        Tuple of (source, subsource, clinical_area).

    Example:
        >>> _parse_source_key("otg__finngen__metabolic")
        ('otg', 'finngen', 'metabolic')
        >>> _parse_source_key("all__all__all")
        ('all', 'all', 'all')
    """
    parts = key.split('__')
    if len(parts) != 3:
        raise ValueError(
            f"Expected key format 'source__subsource__area', got '{key}'"
        )
    return parts[0], parts[1], parts[2]


def write_unified_forest_json(
    forest_results: Dict[str, Dict[str, pd.DataFrame]],
    output_file: str,
    is_baseline: bool = False,
    existing_records: List[Dict] = None,
) -> List[Dict]:
    """
    Write a single structured JSON file with all forest data.

    Output format is a JSON array of objects, each containing:
    - source, subsource, clinical_area, category_type, is_baseline
    - phases: list of per-phase records

    Args:
        forest_results: Dict[source_key, Dict[phase_type, DataFrame]].
        output_file: Path to output JSON file.
        is_baseline: Whether these results are baseline.
        existing_records: Optional list to append to (for combining
            main + baseline into one file).

    Returns:
        The list of records written (for chaining with baseline).
    """
    records: List[Dict] = list(existing_records) if existing_records else []

    for source_key, phase_dict in forest_results.items():
        source, subsource, area = _parse_source_key(source_key)
        for category_type, forest_df in phase_dict.items():
            entry = {
                'source': source,
                'subsource': subsource,
                'clinical_area': area,
                'category_type': category_type,
                'is_baseline': is_baseline,
                'phases': forest_df.to_dict('records'),
            }
            records.append(entry)

    with open(output_file, 'w') as f:
        json.dump(records, f, indent=2)

    logging.info(f"Wrote unified forest JSON: {output_file} ({len(records)} entries)")
    return records


def write_rs_summary_long(
    rs_results: Dict[str, Dict[Tuple[str, str], float]],
    output_file: str,
    is_baseline: bool = False,
    existing_rows: List[Dict] = None,
) -> List[Dict]:
    """
    Write RS summary in long CSV format.

    Columns: source, subsource, clinical_area, category_type, metric, value, is_baseline

    Args:
        rs_results: Dict[source_key, Dict[(phase_type, metric), value]].
        output_file: Path to output CSV file.
        is_baseline: Whether these results are baseline.
        existing_rows: Optional list to append to (for combining
            main + baseline into one file).

    Returns:
        The list of row dicts written (for chaining with baseline).
    """
    rows: List[Dict] = list(existing_rows) if existing_rows else []

    for source_key, rs_dict in rs_results.items():
        source, subsource, area = _parse_source_key(source_key)
        for (category_type, metric), value in rs_dict.items():
            rows.append(
                {
                    'source': source,
                    'subsource': subsource,
                    'clinical_area': area,
                    'category_type': category_type,
                    'metric': metric,
                    'value': value,
                    'is_baseline': is_baseline,
                }
            )

    fieldnames = [
        'source',
        'subsource',
        'clinical_area',
        'category_type',
        'metric',
        'value',
        'is_baseline',
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logging.info(f"Wrote RS summary long CSV: {output_file} ({len(rows)} rows)")
    return rows
