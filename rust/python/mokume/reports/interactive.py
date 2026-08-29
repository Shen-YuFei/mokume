"""
Interactive HTML report for differential expression results.

Generates a standalone HTML file with:
- Interactive volcano plot (hover, click to select proteins)
- Per-protein expression bar chart (updates on selection)
- Searchable/sortable DE results table
- Summary statistics
"""

import json
from string import Template
from typing import Optional

import numpy as np
import pandas as pd

from mokume.core.logger import get_logger

logger = get_logger("mokume.reports.interactive")


def _experimental_samples(protein_df: pd.DataFrame, sample_to_condition: dict) -> list:
    """Return protein column and samples that have a condition mapping."""
    protein_col = protein_df.columns[0]
    sample_cols = [c for c in protein_df.columns if c != protein_col]
    return [
        sample for sample in sample_cols if sample_to_condition.get(sample) is not None
    ]


def _prepare_de_results(de_results: pd.DataFrame) -> pd.DataFrame:
    """Copy DE results and add volcano y-axis values."""
    de = de_results.copy()
    de["neg_log10_p"] = -np.log10(de["adj_pvalue"].clip(1e-300))
    return de


def _condition_colors(sample_to_condition: dict, exp_samples: list) -> dict:
    """Assign stable display colors to conditions."""
    palette = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
    ]
    conditions = sorted(set(sample_to_condition[sample] for sample in exp_samples))
    return {
        condition: palette[index % len(palette)]
        for index, condition in enumerate(conditions)
    }


def _summary_counts(de: pd.DataFrame) -> dict:
    """Return DE significance summary counts."""
    return {
        "n_total": len(de),
        "n_up": (de["significance"] == "UP").sum(),
        "n_down": (de["significance"] == "DOWN").sum(),
        "n_unchanged": (de["significance"] == "Unchanged").sum(),
        "n_not_tested": (de["significance"] == "NotTested").sum(),
    }


def _intensity_payload(protein_df: pd.DataFrame, exp_samples: list) -> dict:
    """Build per-protein intensity payload for the interactive chart."""
    protein_col = protein_df.columns[0]
    intensity_data = {}
    for _, row in protein_df.iterrows():
        intensities = {}
        for sample in exp_samples:
            val = row.get(sample)
            if pd.notna(val) and val > 0:
                intensities[sample] = float(val)
        if intensities:
            intensity_data[str(row[protein_col])] = intensities
    return intensity_data


def _table_rows(de: pd.DataFrame, highlight_genes: list) -> list[dict]:
    """Build DE result rows for the interactive table."""
    return [
        {
            "protein": str(row["ProteinName"]),
            "log2FC": round(float(row["log2FC"]), 4),
            "pvalue": f"{float(row['adj_pvalue']):.2e}",
            "pvalue_raw": float(row["adj_pvalue"]),
            "sig": str(row["significance"]),
            "highlighted": str(row["ProteinName"]) in highlight_genes,
        }
        for _, row in de.iterrows()
    ]


def _volcano_traces(de: pd.DataFrame) -> list[dict]:
    """Build volcano point traces grouped by significance class."""
    traces = []
    for sig_class, color, name in [
        ("Unchanged", "#cccccc", "Unchanged"),
        ("NotTested", "#666666", "NotTested"),
        ("UP", "#d62728", "UP"),
        ("DOWN", "#1f77b4", "DOWN"),
    ]:
        subset = de[de["significance"] == sig_class]
        if subset.empty:
            continue
        traces.append(
            {
                "x": subset["log2FC"].tolist(),
                "y": subset["neg_log10_p"].tolist(),
                "text": subset["ProteinName"].tolist(),
                "customdata": [
                    [
                        row["ProteinName"],
                        f"{row['log2FC']:.3f}",
                        f"{row['adj_pvalue']:.2e}",
                        row["significance"],
                    ]
                    for _, row in subset.iterrows()
                ],
                "color": color,
                "name": f"{name} ({len(subset)})",
            }
        )
    return traces


def _highlight_annotations(de: pd.DataFrame, highlight_genes: list) -> list[dict]:
    """Build annotation payload for highlighted proteins."""
    annotations = []
    for accession in highlight_genes:
        match = de[de["ProteinName"] == accession]
        if match.empty:
            continue
        row = match.iloc[0]
        annotations.append(
            {
                "x": float(row["log2FC"]),
                "y": float(row["neg_log10_p"]),
                "text": str(accession),
            }
        )
    return annotations


def generate_de_report(
    de_results: pd.DataFrame,
    protein_df: pd.DataFrame,
    sample_to_condition: dict,
    output_html: str,
    title: str = "Differential Expression Report",
    highlight_genes: Optional[list] = None,
    log2fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
) -> str:
    """
    Generate an interactive HTML report for differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        DE results with columns: ProteinName, log2FC, adj_pvalue, significance.
    protein_df : pd.DataFrame
        Wide-format protein intensity matrix (ProteinName | sample1 | sample2 | ...).
    sample_to_condition : dict
        Mapping from sample name to condition.
    output_html : str
        Output path for the HTML report.
    title : str
        Report title.
    highlight_genes : list, optional
        Protein accessions to highlight/label on the volcano plot.
    log2fc_threshold : float
        log2FC threshold used for significance.
    fdr_threshold : float
        FDR threshold used for significance.

    Returns
    -------
    str
        Path to the generated HTML file.
    """
    if highlight_genes is None:
        highlight_genes = []

    exp_samples = _experimental_samples(protein_df, sample_to_condition)
    de = _prepare_de_results(de_results)
    counts = _summary_counts(de)

    html = _build_html(
        title=title,
        volcano_traces=_volcano_traces(de),
        annotations=_highlight_annotations(de, highlight_genes),
        intensity_data=json.dumps(_intensity_payload(protein_df, exp_samples)),
        table_rows=json.dumps(_table_rows(de, highlight_genes)),
        sample_to_condition=json.dumps(
            {s: sample_to_condition[s] for s in exp_samples}
        ),
        cond_colors=json.dumps(_condition_colors(sample_to_condition, exp_samples)),
        exp_samples=json.dumps(exp_samples),
        n_total=counts["n_total"],
        n_up=counts["n_up"],
        n_down=counts["n_down"],
        n_unchanged=counts["n_unchanged"],
        n_not_tested=counts["n_not_tested"],
        log2fc_threshold=log2fc_threshold,
        fdr_threshold=fdr_threshold,
    )

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info("Interactive report saved to %s", output_html)
    return output_html


def _build_html(
    title,
    volcano_traces,
    annotations,
    intensity_data,
    table_rows,
    sample_to_condition,
    cond_colors,
    exp_samples,
    n_total,
    n_up,
    n_down,
    n_unchanged,
    n_not_tested,
    log2fc_threshold,
    fdr_threshold,
):
    """Build the complete HTML report string."""
    # Build plotly volcano trace JSON
    plotly_traces = []
    for t in volcano_traces:
        plotly_traces.append(
            {
                "type": "scatter",
                "mode": "markers",
                "x": t["x"],
                "y": t["y"],
                "text": t["text"],
                "customdata": t["customdata"],
                "marker": {"color": t["color"], "size": 5, "opacity": 0.7},
                "name": t["name"],
                "hovertemplate": (
                    "<b>%{customdata[0]}</b><br>"
                    "log2FC: %{customdata[1]}<br>"
                    "adj p-value: %{customdata[2]}<br>"
                    "Significance: %{customdata[3]}<extra></extra>"
                ),
            }
        )

    volcano_layout = {
        "title": {"text": f"Volcano Plot: {title}", "font": {"size": 16}},
        "xaxis": {"title": "log2 Fold Change", "zeroline": True},
        "yaxis": {"title": "-log10(adjusted p-value)"},
        "hovermode": "closest",
        "legend": {"x": 0.01, "y": 0.99},
        "shapes": [
            {
                "type": "line",
                "x0": log2fc_threshold,
                "x1": log2fc_threshold,
                "y0": 0,
                "y1": 1,
                "yref": "paper",
                "line": {"dash": "dash", "color": "grey", "width": 1},
            },
            {
                "type": "line",
                "x0": -log2fc_threshold,
                "x1": -log2fc_threshold,
                "y0": 0,
                "y1": 1,
                "yref": "paper",
                "line": {"dash": "dash", "color": "grey", "width": 1},
            },
            {
                "type": "line",
                "y0": -np.log10(fdr_threshold),
                "y1": -np.log10(fdr_threshold),
                "x0": 0,
                "x1": 1,
                "xref": "paper",
                "line": {"dash": "dash", "color": "grey", "width": 1},
            },
        ],
        "annotations": [
            {
                "x": a["x"],
                "y": a["y"],
                "text": a["text"],
                "showarrow": True,
                "arrowhead": 2,
                "arrowsize": 0.8,
                "ax": 20,
                "ay": -25,
                "font": {"size": 10, "color": "darkgreen"},
            }
            for a in annotations
        ],
    }

    plotly_traces_json = json.dumps(plotly_traces)
    volcano_layout_json = json.dumps(volcano_layout)

    return Template("""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>$title</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
               background: #f5f5f5; color: #333; }
        .header { background: #2c3e50; color: white; padding: 20px 30px; }
        .header h1 { font-size: 24px; margin-bottom: 8px; }
        .stats { display: flex; gap: 20px; margin-top: 10px; }
        .stat-box { background: rgba(255,255,255,0.15); padding: 8px 16px;
                     border-radius: 6px; text-align: center; }
        .stat-box .number { font-size: 22px; font-weight: bold; }
        .stat-box .label { font-size: 11px; opacity: 0.8; }
        .container { max-width: 1400px; margin: 0 auto; padding: 20px; }
        .row { display: flex; gap: 20px; margin-bottom: 20px; }
        .col-8 { flex: 2; }
        .col-4 { flex: 1; }
        .card { background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                 padding: 16px; }
        .card h3 { margin-bottom: 12px; font-size: 15px; color: #2c3e50; }
        #volcano-plot { width: 100%; height: 500px; }
        #expression-plot { width: 100%; height: 350px; }
        .table-controls { margin-bottom: 10px; display: flex; gap: 10px; align-items: center; }
        .table-controls input { padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;
                                  font-size: 13px; width: 250px; }
        .table-controls select { padding: 6px 8px; border: 1px solid #ddd; border-radius: 4px;
                                   font-size: 13px; }
        #de-table { width: 100%; border-collapse: collapse; font-size: 12px; }
        #de-table th { background: #34495e; color: white; padding: 8px 10px; cursor: pointer;
                        text-align: left; position: sticky; top: 0; }
        #de-table th:hover { background: #2c3e50; }
        #de-table td { padding: 6px 10px; border-bottom: 1px solid #eee; }
        #de-table tr:hover { background: #e8f4fd; cursor: pointer; }
        #de-table tr.selected { background: #d4edda; }
        .sig-up { color: #d62728; font-weight: bold; }
        .sig-down { color: #1f77b4; font-weight: bold; }
        .sig-unchanged { color: #999; }
        .sig-nottested { color: #666; }
        .table-wrapper { max-height: 400px; overflow-y: auto; border: 1px solid #ddd;
                          border-radius: 4px; }
        .info-text { font-size: 12px; color: #666; margin-top: 8px; }
        .highlight { background: #fff3cd; }
    </style>
</head>
<body>

<div class="header">
    <h1>$title</h1>
    <div class="stats">
        <div class="stat-box">
            <div class="number">$n_total</div>
            <div class="label">Proteins Tested</div>
        </div>
        <div class="stat-box" style="border-left: 3px solid #d62728;">
            <div class="number" style="color: #ff6b6b;">$n_up</div>
            <div class="label">Upregulated</div>
        </div>
        <div class="stat-box" style="border-left: 3px solid #1f77b4;">
            <div class="number" style="color: #74b9ff;">$n_down</div>
            <div class="label">Downregulated</div>
        </div>
        <div class="stat-box">
            <div class="number">$n_unchanged</div>
            <div class="label">Unchanged</div>
        </div>
        <div class="stat-box">
            <div class="number">$n_not_tested</div>
            <div class="label">Not tested</div>
        </div>
        <div class="stat-box">
            <div class="number">$log2fc_threshold</div>
            <div class="label">|log2FC| cutoff</div>
        </div>
        <div class="stat-box">
            <div class="number">$fdr_threshold</div>
            <div class="label">FDR cutoff</div>
        </div>
    </div>
</div>

<div class="container">
    <div class="row">
        <div class="col-8">
            <div class="card">
                <h3>Interactive Volcano Plot</h3>
                <div id="volcano-plot"></div>
                <p class="info-text">Click a point to view per-sample expression below.</p>
            </div>
        </div>
        <div class="col-4">
            <div class="card">
                <h3>Protein Expression</h3>
                <div id="expression-plot"></div>
                <p class="info-text" id="expression-info">Select a protein from the volcano plot or table.</p>
            </div>
        </div>
    </div>
    <div class="row">
        <div style="flex: 1;">
            <div class="card">
                <h3>Differential Expression Results</h3>
                <div class="table-controls">
                    <input type="text" id="search-input" placeholder="Search protein accession...">
                    <select id="sig-filter">
                        <option value="all">All</option>
                        <option value="UP">UP only</option>
                        <option value="DOWN">DOWN only</option>
                        <option value="NotTested">Not tested</option>
                        <option value="significant">Significant (UP + DOWN)</option>
                    </select>
                    <span id="table-count" style="font-size: 12px; color: #666;"></span>
                </div>
                <div class="table-wrapper">
                    <table id="de-table">
                        <thead>
                            <tr>
                                <th onclick="sortTable(0)">Protein</th>
                                <th onclick="sortTable(1)">log2FC</th>
                                <th onclick="sortTable(2)">adj p-value</th>
                                <th onclick="sortTable(3)">Significance</th>
                            </tr>
                        </thead>
                        <tbody id="de-tbody"></tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
</div>

<script>
// Data
const intensityData = $intensity_data;
const tableRows = $table_rows;
const sampleToCondition = $sample_to_condition;
const condColors = $cond_colors;
const expSamples = $exp_samples;

// Volcano plot
const volcanoTraces = $plotly_traces_json;
const volcanoLayout = $volcano_layout_json;
volcanoLayout.margin = {t: 50, b: 50, l: 60, r: 20};

Plotly.newPlot('volcano-plot', volcanoTraces, volcanoLayout, {responsive: true});

// Click handler for volcano
document.getElementById('volcano-plot').on('plotly_click', function(data) {
    const protein = data.points[0].customdata[0];
    showExpression(protein);
    highlightTableRow(protein);
});

// Expression plot
function showExpression(protein) {
    const data = intensityData[protein];
    if (!data) {
        document.getElementById('expression-info').textContent = 'No intensity data for ' + protein;
        return;
    }

    // Find DE info
    const deRow = tableRows.find(r => r.protein === protein);
    const fcText = deRow ? ' (log2FC=' + deRow.log2FC.toFixed(2) + ', ' + deRow.sig + ')' : '';

    const samples = expSamples.filter(s => s in data);
    const x = samples;
    const y = samples.map(s => Math.log2(data[s]));
    const colors = samples.map(s => condColors[sampleToCondition[s]] || '#999');

    const trace = {
        type: 'bar',
        x: x,
        y: y,
        marker: { color: colors },
        hovertemplate: '<b>%{x}</b><br>log2 intensity: %{y:.2f}<extra></extra>',
    };

    const layout = {
        title: { text: protein + fcText, font: { size: 13 } },
        xaxis: { tickangle: 45, tickfont: { size: 9 } },
        yaxis: { title: 'log2 Intensity' },
        margin: { t: 40, b: 80, l: 50, r: 10 },
        showlegend: false,
    };

    Plotly.newPlot('expression-plot', [trace], layout, {responsive: true});
    document.getElementById('expression-info').textContent =
        'Showing: ' + protein + ' | Condition colors match legend';
}

// Populate table
function populateTable(filter, search) {
    const tbody = document.getElementById('de-tbody');
    tbody.innerHTML = '';
    let count = 0;

    tableRows.forEach(row => {
        if (filter === 'UP' && row.sig !== 'UP') return;
        if (filter === 'DOWN' && row.sig !== 'DOWN') return;
        if (filter === 'NotTested' && row.sig !== 'NotTested') return;
        if (filter === 'significant' && row.sig !== 'UP' && row.sig !== 'DOWN') return;
        if (search && !row.protein.toLowerCase().includes(search.toLowerCase())) return;

        count++;
        const tr = document.createElement('tr');
        tr.id = 'row-' + row.protein;
        tr.onclick = function() {
            showExpression(row.protein);
            highlightTableRow(row.protein);
        };

        const sigClass = row.sig === 'UP' ? 'sig-up'
            : row.sig === 'DOWN' ? 'sig-down'
            : row.sig === 'NotTested' ? 'sig-nottested'
            : 'sig-unchanged';

        tr.innerHTML = '<td>' + (row.highlighted ? '<b>' + row.protein + '</b>' : row.protein) + '</td>'
            + '<td style="text-align:right">' + row.log2FC.toFixed(3) + '</td>'
            + '<td style="text-align:right">' + row.pvalue + '</td>'
            + '<td class="' + sigClass + '">' + row.sig + '</td>';
        tbody.appendChild(tr);
    });

    document.getElementById('table-count').textContent = count + ' / ' + tableRows.length + ' proteins';
}

function highlightTableRow(protein) {
    document.querySelectorAll('#de-table tr.selected').forEach(tr => tr.classList.remove('selected'));
    const row = document.getElementById('row-' + protein);
    if (row) {
        row.classList.add('selected');
        row.scrollIntoView({ behavior: 'smooth', block: 'center' });
    }
}

// Sort
let sortDir = [1, 1, 1, 1];
function sortTable(col) {
    sortDir[col] *= -1;
    tableRows.sort((a, b) => {
        const keys = ['protein', 'log2FC', 'pvalue_raw', 'sig'];
        let va = a[keys[col]], vb = b[keys[col]];
        if (typeof va === 'string') return va.localeCompare(vb) * sortDir[col];
        return (va - vb) * sortDir[col];
    });
    const filter = document.getElementById('sig-filter').value;
    const search = document.getElementById('search-input').value;
    populateTable(filter, search);
}

// Filters
document.getElementById('search-input').addEventListener('input', function() {
    const filter = document.getElementById('sig-filter').value;
    populateTable(filter, this.value);
});
document.getElementById('sig-filter').addEventListener('change', function() {
    const search = document.getElementById('search-input').value;
    populateTable(this.value, search);
});

// Initial render
populateTable('all', '');

// Show first highlighted gene if available
const firstHighlight = tableRows.find(r => r.highlighted);
if (firstHighlight) showExpression(firstHighlight.protein);
</script>

<div style="text-align: center; padding: 20px; color: #999; font-size: 11px;">
    Generated by <a href="https://github.com/bigbio/mokume" style="color: #666;">mokume</a>
</div>
</body>
</html>""").substitute(
        title=title,
        n_total=n_total,
        n_up=n_up,
        n_down=n_down,
        n_unchanged=n_unchanged,
        n_not_tested=n_not_tested,
        log2fc_threshold=log2fc_threshold,
        fdr_threshold=fdr_threshold,
        intensity_data=intensity_data,
        table_rows=table_rows,
        sample_to_condition=sample_to_condition,
        cond_colors=cond_colors,
        exp_samples=exp_samples,
        plotly_traces_json=plotly_traces_json,
        volcano_layout_json=volcano_layout_json,
    )
