import os
import json
import pandas as pd
from string import Template
from pathlib import Path

def write_plot(template_name, output_file, **template_vars):
    """
    Write an interactive HTML plot using a template file.

    This function loads an HTML template from the templates directory,
    fills in variables using Python's string.Template, and writes the
    result to an output file.

    Args:
        template_name: Name of the template file (e.g., 'diversity_heatmaps.html')
        output_file: Path to output HTML file
        **template_vars: Keyword arguments for template variables

    Example:
        write_plot(
            'diversity_heatmaps.html',
            'output.html',
            title='My Analysis',
            stats_json=json.dumps(data),
            pair_options='<option>A</option>'
        )
    """
    # Get the template directory (relative to this file)
    templates_dir = Path(__file__).parent / 'templates'
    template_path = templates_dir / template_name

    # Check if template exists
    if not template_path.exists():
        raise FileNotFoundError(f"Template not found: {template_path}")

    # Read the template
    with open(template_path, 'r') as f:
        template_content = f.read()

    # Create Template object and substitute variables
    template = Template(template_content)
    html_output = template.safe_substitute(**template_vars)

    # Write the output
    with open(output_file, 'w') as f:
        f.write(html_output)


def save_figures(figs, output_file, title="Interactive Plots", table_of_contents=False):
    """ A nice function to save a list of figures to a single html file """
    html_content = f"""
    <html>
    <head>
        <title>{title}</title>
    </head>
    <body>
        <div style="margin-bottom: 20px;">
            {figs[0].to_html(include_plotlyjs="cdn", full_html=False)}
        </div>
    """
    for fig in figs[1:]:
        html_content += f"""
        <div style="margin-bottom: 20px;">
            {fig.to_html(include_plotlyjs=False, full_html=False)}
        </div>
        """
    html_content += """
    </body>
    </html>
    """
    with open(output_file, 'w') as f:
        f.write(html_content)

def make_table_of_figures(
    figs,
    output_file,
    title="Interactive Plots",
    title_delimiter=None,
    column_titles=None,
    ):
    """ Make a table of figures in the LAP UI, by saving figure in different files, but returning links to each figure in a table
    Args:
        figs: List of figures
        output_file: Output file of the table
        title: Title of the table
    """
    filepath_no_extension = os.path.splitext(output_file)[0]
    base_filename = os.path.basename(filepath_no_extension)
    html_content = f"""
    <html>
    <head>
        <title>{title}</title>
        <style>
            table {{ border-collapse: collapse; width: 80%; margin: 20px auto; }}
            th, td {{ border: 1px solid #888; padding: 8px 12px; text-align: left; }}
            th {{ background-color: #f2f2f2; cursor: pointer; }}
            tr:nth-child(even) {{ background-color: #fafafa; }}
            h1 {{ text-align: center; }}
            #filterInput {{ margin: 20px auto; display: block; width: 60%; padding: 8px; font-size: 1em; }}
        </style>
        <script>
        // Filtering
        function filterTable() {{
            var input = document.getElementById('filterInput');
            var filter = input.value.toLowerCase();
            var table = document.getElementById('figTable');
            var trs = table.getElementsByTagName('tr');
            for (var i = 1; i < trs.length; i++) {{
                var tds = trs[i].getElementsByTagName('td');
                var show = false;
                for (var j = 0; j < tds.length; j++) {{
                    if (tds[j].textContent.toLowerCase().indexOf(filter) > -1) {{
                        show = true;
                        break;
                    }}
                }}
                trs[i].style.display = show ? '' : 'none';
            }}
        }}
        // Sorting
        function sortTable(n) {{
            var table = document.getElementById('figTable');
            var switching = true;
            var dir = 'asc';
            var switchcount = 0;
            while (switching) {{
                switching = false;
                var rows = table.rows;
                for (var i = 1; i < (rows.length - 1); i++) {{
                    var shouldSwitch = false;
                    var x = rows[i].getElementsByTagName('TD')[n];
                    var y = rows[i + 1].getElementsByTagName('TD')[n];
                    if (dir == 'asc') {{
                        if (x.textContent.toLowerCase() > y.textContent.toLowerCase()) {{
                            shouldSwitch = true;
                            break;
                        }}
                    }} else if (dir == 'desc') {{
                        if (x.textContent.toLowerCase() < y.textContent.toLowerCase()) {{
                            shouldSwitch = true;
                            break;
                        }}
                    }}
                }}
                if (shouldSwitch) {{
                    rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                    switching = true;
                    switchcount++;
                }} else {{
                    if (switchcount == 0 && dir == 'asc') {{
                        dir = 'desc';
                        switching = true;
                    }}
                }}
            }}
        }}
        </script>
    </head>
    <body>
        <h1>{title}</h1>
        <input type="text" id="filterInput" onkeyup="filterTable()" placeholder="Filter table..." />
        <table id="figTable">
    """
    # Construct table headers
    if column_titles is None:
        column_titles = ['Figure Title', 'Link']
        custom_title = False
    else:
        column_titles = column_titles + ['Link']
        custom_title = True
    header_html = ""
    for idx, title in enumerate(column_titles):
        header_html += f"""
            <th onclick='sortTable({idx})'>{title}</th>
        """
    html_content += f"""
        <tr>
            {header_html}
        </tr>
    """
    for i, fig in enumerate(figs):
        fig_file = f"{filepath_no_extension}_{i}.html"
        fig_url = f"{base_filename}_{i}.html"
        fig.write_html(fig_file)
        if custom_title:
            title_parts = fig.layout.title.text.split(title_delimiter)
        else:
            title_parts = [fig.layout.title.text]
        row_html = ""
        for title in title_parts:
            row_html += f"""
                <td>{title}</td>
            """
        row_html += f"""
            <td><a href='{fig_url}'>link</a></td>
        """
        html_content += f"""
            <tr>
                {row_html}
            </tr>
        """
    html_content += """
        </table>
    </body>
    </html>
    """
    with open(output_file, 'w') as f:
        f.write(html_content)


def create_interactive_diversity_figure(stats_df, output_file, title="Target Diversity Analysis"):
    """
    Create a self-contained interactive HTML file for target diversity analysis.

    This function embeds all data as JSON and uses JavaScript to handle state management
    and dynamic Plotly graph rendering. Everything is contained in a single HTML file
    suitable for static hosting in the LAP UI.

    Args:
        stats_df: DataFrame with columns: pair, run, other, category_type, area, phase,
                  analysis, analysis_type, count, percent
        output_file: Path to output HTML file
        title: Title for the page
    """
    # Convert dataframe to JSON
    stats_json = stats_df.to_json(orient='records')

    # Get unique values for dropdowns
    runs = sorted(stats_df['run'].unique().tolist())
    phases = ['all', 'Preclinical', 'Phase I', 'Phase II', 'Phase III', 'Launched']
    areas = sorted([area for area in stats_df['area'].unique() if area != 'all'])
    all_areas = ['all'] + areas

    # Generate dropdown options HTML
    run_options = '\n'.join(f'<option value="{run}">{run}</option>' for run in runs)

    # Use the templating system
    write_plot(
        'diversity_heatmaps.html',
        output_file,
        title=title,
        stats_json=stats_json,
        phase_order_json=json.dumps(phases),
        all_areas_json=json.dumps(all_areas),
        run_options=run_options
    )