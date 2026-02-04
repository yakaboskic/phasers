#!/usr/bin/env python3
"""
Phasers CLI - Command-line interface for relative success calculations.

Usage:
    phasers continuous --help
    phasers plot --help
    phasers summary --help
"""

import argparse
import logging
import sys
import traceback

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from tqdm import tqdm

tqdm.pandas()

from .core import (
    calculate_continuous_rs_incremental,
    write_continuous_rs_points,
    load_and_merge_all_data,
    calculate_all_similarities,
    write_merged_data,
)
from .plotting import save_figures


def handle_continuous(args):
    """
    Calculate continuous RS curve using ROC-style incremental updates.

    This analyzes how RS changes as you select T-I pairs by score
    using efficient incremental count updates.
    """
    logging.info("=" * 60)
    logging.info("PHASERS: CONTINUOUS RS CALCULATION")
    logging.info("=" * 60)

    # Load data
    merged_df = None
    if hasattr(args, 'merged_data_file') and args.merged_data_file:
        logging.info(f"Loading pre-merged data from: {args.merged_data_file}")
        compression = 'gzip' if args.merged_data_file.endswith('.gz') else None
        merged_df = pd.read_csv(args.merged_data_file, sep='\t', compression=compression, low_memory=False)
        logging.info(f"Loaded {len(merged_df)} rows from pre-merged file")

    # Parse score columns (comma-separated)
    # Format: column_name or column_name:ascending
    score_columns = []
    score_ascending_map = {}
    for col_spec in args.score_column.split(','):
        col_spec = col_spec.strip()
        if ':' in col_spec:
            col_name, direction = col_spec.rsplit(':', 1)
            col_name = col_name.strip()
            if direction.lower() == 'ascending':
                score_ascending_map[col_name] = True
                logging.info(f"Score column '{col_name}' will be sorted ascending (lower is better)")
            else:
                logging.warning(f"Unknown direction '{direction}' for column '{col_name}', using default descending")
                score_ascending_map[col_name] = False
        else:
            col_name = col_spec
            score_ascending_map[col_name] = False
        score_columns.append(col_name)
    logging.info(f"Requested score columns: {score_columns}")

    # Load data if not pre-merged
    if merged_df is None:
        logging.info("Loading data from individual files...")
        merged_df, _ = load_and_merge_all_data(
            pharmaprojects_file=getattr(args, 'pharmaprojects_file', None),
            indications_file=getattr(args, 'indications_file', None),
            associations_file=getattr(args, 'associations_file', None),
            similarity_matrix_file=getattr(args, 'similarity_matrix_file', None),
            similarity_matrix_format=args.similarity_matrix_format,
            no_genetic_insight=args.no_genetic_insight,
            association_thresholds=None,
        )

    # Filter to source to check available columns
    merged_df['source'] = merged_df['source'].str.lower()
    source_data = merged_df[merged_df['source'] == args.source_name.lower()]
    available_columns = set(source_data.columns)

    # Apply filters if specified
    if args.filter:
        from .core import apply_association_filters_to_merged_data
        logging.info(f"Applying filters: {args.filter}")
        merged_df = apply_association_filters_to_merged_data(merged_df, args.filter)
        logging.info(f"Filtered data to {len(merged_df)} rows")

    # Match requested columns to available columns
    matched_columns = [col for col in score_columns if col in available_columns]
    missing_columns = [col for col in score_columns if col not in available_columns]

    if missing_columns:
        logging.warning(f"Skipping missing columns: {missing_columns}")

    if not matched_columns:
        logging.error(f"No valid score columns found! Requested: {score_columns}, Available: {list(available_columns)}")
        sys.exit(1)

    logging.info(f"Running continuous RS for {len(matched_columns)} score columns: {matched_columns}")

    # Calculate continuous RS for each score column
    all_results = []
    for score_col in matched_columns:
        logging.info(f"\n{'='*60}")
        logging.info(f"Processing score column: {score_col}")
        logging.info(f"{'='*60}")

        try:
            is_ascending = score_ascending_map.get(score_col, False)

            results_df = calculate_continuous_rs_incremental(
                score_column=score_col,
                source_name=args.source_name,
                phase_type=args.phase_type,
                similarity_threshold=args.similarity_matrix_threshold,
                score_strategy=args.score_strategy,
                score_ascending=is_ascending,
                merged_df=merged_df,
                pharmaprojects_file=None,
                indications_file=None,
                associations_file=None,
                similarity_matrix_file=None,
                similarity_matrix_format=args.similarity_matrix_format,
                no_genetic_insight=args.no_genetic_insight,
                deduplicated_data_file=None,
            )

            results_df['score_column'] = score_col
            all_results.append(results_df)
            logging.info(f"Successfully calculated RS for {score_col}: {len(results_df)} points")

        except Exception as e:
            logging.error(f"Failed to calculate RS for column '{score_col}': {e}")
            logging.error(traceback.format_exc())
            continue

    if not all_results:
        logging.error("No score columns processed successfully!")
        sys.exit(1)

    # Combine and write results
    combined_results = pd.concat(all_results, ignore_index=True)
    write_continuous_rs_points(combined_results, args.output_file)

    logging.info("=" * 60)
    logging.info("CONTINUOUS RS CALCULATION COMPLETE")
    logging.info(f"Output: {args.output_file}")
    logging.info("=" * 60)


def handle_plot(args):
    """
    Plot continuous RS curve from points file.
    """
    logging.info("=" * 60)
    logging.info("PHASERS: PLOTTING RS CURVE")
    logging.info("=" * 60)

    # Load data
    df = pd.read_csv(args.input_file)
    logging.info(f"Loaded {len(df)} points from {args.input_file}")

    # Filter by alpha range
    df = df[(df['alpha'] >= args.start_alpha) & (df['alpha'] <= args.max_alpha)]
    logging.info(f"Filtered to alpha range [{args.start_alpha}, {args.max_alpha}]: {len(df)} points")

    # Filter by score columns if specified
    if args.score_column_filter:
        requested_cols = [col.strip() for col in args.score_column_filter.split(',')]
        df = df[df['score_column'].isin(requested_cols)].copy()
        logging.info(f"Filtered to score columns {requested_cols}: {len(df)} points")

    # Create figure
    fig = go.Figure()

    # Color palette for score columns
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    score_columns = df['score_column'].unique() if 'score_column' in df.columns else ['default']

    for idx, score_col in enumerate(score_columns):
        if 'score_column' in df.columns:
            df_score = df[df['score_column'] == score_col].sort_values('alpha')
        else:
            df_score = df.sort_values('alpha')

        color = colors[idx % len(colors)]

        # Add CI band if available
        if 'RS_cum_lower_ci' in df_score.columns and 'RS_cum_upper_ci' in df_score.columns:
            fig.add_trace(go.Scatter(
                x=df_score['alpha'],
                y=df_score['RS_cum_upper_ci'],
                mode='lines',
                line=dict(width=0),
                showlegend=False,
                hoverinfo='skip',
                legendgroup=score_col
            ))
            fig.add_trace(go.Scatter(
                x=df_score['alpha'],
                y=df_score['RS_cum_lower_ci'],
                mode='lines',
                line=dict(width=0),
                fill='tonexty',
                fillcolor=f'rgba{tuple(list(int(color[i:i+2], 16) for i in (1, 3, 5)) + [0.2])}',
                showlegend=False,
                hoverinfo='skip',
                legendgroup=score_col
            ))

        # Main RS curve
        fig.add_trace(go.Scatter(
            x=df_score['alpha'],
            y=df_score['RS_cum'],
            mode='lines',
            name=score_col,
            line=dict(color=color, width=2),
            hovertemplate=f'<b>{score_col}</b><br>Alpha: %{{x:.4f}}<br>RS_cum: %{{y:.4f}}<extra></extra>',
            legendgroup=score_col
        ))

    # Add baseline
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray",
                  annotation_text="RS = 1.0 (baseline)", annotation_position="right")

    # Update layout
    fig.update_layout(
        title="Continuous Relative Success Curve",
        xaxis_title="Alpha (Fraction of T-I Pairs Selected)",
        yaxis_title="Cumulative Relative Success (RS_cum)",
        xaxis_range=[args.start_alpha, args.max_alpha],
        hovermode='x unified',
        template='plotly_white',
        showlegend=True,
    )

    # Save
    save_figures([fig], args.output_file, title="Phasers RS Curve")
    logging.info(f"Saved plot to: {args.output_file}")


def handle_summary(args):
    """
    Generate summary statistics from RS points file.
    """
    logging.info("=" * 60)
    logging.info("PHASERS: GENERATING SUMMARY STATISTICS")
    logging.info("=" * 60)

    # Load data
    df = pd.read_csv(args.input_file)
    logging.info(f"Loaded {len(df)} points from {args.input_file}")

    # Extract summary metrics (they're the same for all rows within a score_column)
    summary_cols = ['EA_RS', 'MLRS', 'RS_max', 'RS_at_10pct']
    available_summary_cols = [col for col in summary_cols if col in df.columns]

    if 'score_column' in df.columns:
        summary_df = df.groupby('score_column')[available_summary_cols].first().reset_index()
    else:
        summary_df = df[available_summary_cols].head(1)

    # Add source info if available
    if 'source' in df.columns:
        summary_df['source'] = df['source'].iloc[0]

    # Write summary
    summary_df.to_csv(args.output_file, index=False)
    logging.info(f"Saved summary to: {args.output_file}")

    # Print summary
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))


def handle_merge(args):
    """
    Merge input data files and calculate similarities.

    This creates a merged data file that can be reused for multiple
    RS calculations, significantly speeding up analysis.
    """
    logging.info("=" * 60)
    logging.info("PHASERS: MERGING DATA FILES")
    logging.info("=" * 60)

    # Load and merge all data
    logging.info("Loading and merging data files...")
    merged_df, similarity_matrix = load_and_merge_all_data(
        pharmaprojects_file=args.pharmaprojects_file,
        indications_file=args.indications_file,
        associations_file=args.associations_file,
        similarity_matrix_file=args.similarity_matrix_file,
        pigean_file=getattr(args, 'pigean_file', None),
        similarity_matrix_format=args.similarity_matrix_format,
        no_genetic_insight=args.no_genetic_insight,
        association_thresholds=args.filter if args.filter else None,
    )

    logging.info(f"Merged data: {len(merged_df)} rows")

    # Calculate similarities
    logging.info("Calculating similarities...")
    merged_df = calculate_all_similarities(merged_df, similarity_matrix)

    # Save merged data
    logging.info(f"Saving merged data to: {args.output_file}")
    compression = 'gzip' if args.output_file.endswith('.gz') else None
    merged_df.to_csv(args.output_file, sep='\t', index=False, compression=compression)

    logging.info("=" * 60)
    logging.info("MERGE COMPLETE")
    logging.info(f"Output: {args.output_file}")
    logging.info(f"Total rows: {len(merged_df)}")
    logging.info(f"Sources: {merged_df['source'].dropna().unique().tolist()}")
    logging.info("=" * 60)


def main():
    """Main CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s] %(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    parser = argparse.ArgumentParser(
        prog='phasers',
        description='Phasers - Clinical Phase Relative Success Calculator'
    )
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')

    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands')

    # === CONTINUOUS COMMAND ===
    parser_continuous = subparsers.add_parser(
        'continuous',
        help='Calculate continuous RS curve using incremental updates'
    )

    # Input data
    input_group = parser_continuous.add_argument_group('input data')
    input_group.add_argument(
        '--merged-data-file', type=str,
        help='Pre-merged data file (TSV or TSV.gz) - FAST option'
    )
    input_group.add_argument(
        '--pharmaprojects-file', type=str,
        help='Pharmaprojects file (TSV) - required if not using --merged-data-file'
    )
    input_group.add_argument(
        '--indications-file', type=str,
        help='Indications file (TSV) - required if not using --merged-data-file'
    )
    input_group.add_argument(
        '--associations-file', type=str,
        help='Associations file (TSV) - required if not using --merged-data-file'
    )
    input_group.add_argument(
        '--similarity-matrix-file', type=str,
        help='Similarity matrix file (CSV) - required if not using --merged-data-file'
    )

    # Required
    parser_continuous.add_argument(
        '--output-file', '-o', required=True,
        help='Output CSV file with RS curve points'
    )
    parser_continuous.add_argument(
        '--score-column', '-s', required=True,
        help='Comma-separated score columns. Use column:ascending for metrics where lower is better (e.g., "combined,indirect,p_value:ascending")'
    )
    parser_continuous.add_argument(
        '--source-name', required=True,
        help='Association source name to filter (e.g., "pigean", "falcon")'
    )

    # Optional
    parser_continuous.add_argument(
        '--filter', type=str, default=None, action='append',
        help='Filter to apply (e.g., "falcon:p_value:<=:1e-8:")'
    )
    parser_continuous.add_argument(
        '--phase-type', type=str, default='ccat', choices=['ccat', 'hcat', 'acat'],
        help='Phase type to analyze (default: ccat)'
    )
    parser_continuous.add_argument(
        '--similarity-matrix-format', type=str, default='long', choices=['long', 'matrix'],
        help='Similarity matrix format (default: long)'
    )
    parser_continuous.add_argument(
        '--similarity-matrix-threshold', type=float, default=0.8,
        help='Similarity threshold for genetic support (default: 0.8)'
    )
    parser_continuous.add_argument(
        '--score-strategy', type=str, default='avg', choices=['max', 'min', 'avg'],
        help='Strategy for handling duplicate scores (default: avg)'
    )
    parser_continuous.add_argument(
        '--no-genetic-insight', action='store_true',
        help='Disable genetic insight filter'
    )

    # === PLOT COMMAND ===
    parser_plot = subparsers.add_parser(
        'plot',
        help='Plot continuous RS curve from points file'
    )
    parser_plot.add_argument(
        '--input-file', '-i', required=True,
        help='Input CSV file with RS points (from continuous command)'
    )
    parser_plot.add_argument(
        '--output-file', '-o', required=True,
        help='Output HTML file with interactive plot'
    )
    parser_plot.add_argument(
        '--start-alpha', type=float, default=0.0,
        help='Minimum alpha to plot (default: 0.0)'
    )
    parser_plot.add_argument(
        '--max-alpha', type=float, default=1.0,
        help='Maximum alpha to plot (default: 1.0)'
    )
    parser_plot.add_argument(
        '--score-column-filter', type=str, default=None,
        help='Filter to specific score columns (comma-separated)'
    )

    # === SUMMARY COMMAND ===
    parser_summary = subparsers.add_parser(
        'summary',
        help='Generate summary statistics from RS points file'
    )
    parser_summary.add_argument(
        '--input-file', '-i', required=True,
        help='Input CSV file with RS points'
    )
    parser_summary.add_argument(
        '--output-file', '-o', required=True,
        help='Output CSV file with summary statistics'
    )

    # === MERGE COMMAND ===
    parser_merge = subparsers.add_parser(
        'merge',
        help='Merge input data files and calculate similarities (creates reusable merged data file)'
    )
    parser_merge.add_argument(
        '--pharmaprojects-file', required=True,
        help='Pharmaprojects file (TSV) - drug-target-indication clinical trial data'
    )
    parser_merge.add_argument(
        '--indications-file', required=True,
        help='Indications file (TSV) - indication metadata with genetic insight'
    )
    parser_merge.add_argument(
        '--associations-file', required=True,
        help='Associations file (TSV) - gene-indication association scores'
    )
    parser_merge.add_argument(
        '--similarity-matrix-file', required=True,
        help='Similarity matrix file (CSV) - MeSH ID similarity scores'
    )
    parser_merge.add_argument(
        '--output-file', '-o', required=True,
        help='Output merged data file (TSV or TSV.gz)'
    )
    parser_merge.add_argument(
        '--pigean-file', type=str, default=None,
        help='Optional separate PIGEAN associations file'
    )
    parser_merge.add_argument(
        '--similarity-matrix-format', type=str, default='long', choices=['long', 'matrix'],
        help='Similarity matrix format (default: long)'
    )
    parser_merge.add_argument(
        '--no-genetic-insight', action='store_true',
        help='Disable genetic insight filter (include all indications)'
    )
    parser_merge.add_argument(
        '--filter', type=str, default=None, action='append',
        help='Filter to apply to associations (e.g., "falcon:p_value:<=:1e-8:")'
    )

    # Parse and route
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        if args.command == 'merge':
            handle_merge(args)
        elif args.command == 'continuous':
            handle_continuous(args)
        elif args.command == 'plot':
            handle_plot(args)
        elif args.command == 'summary':
            handle_summary(args)
        else:
            parser.print_help()
            sys.exit(1)

    except KeyboardInterrupt:
        logging.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
