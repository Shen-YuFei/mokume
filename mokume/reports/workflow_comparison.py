"""
Multi-workflow comparison report for proteomics quantification.

Compares quantification workflows side-by-side using:
- PCA plots per workflow
- Summary metrics table
- DE concordance (overlap, log2FC correlation, concordance-at-top)
- Quality metrics radar chart
- Known marker gene panel

Designed to answer: "Which workflow preserves biology better?"
"""

import html
from typing import Optional

import numpy as np

from mokume.reports.qc_report import compute_qc_metrics, _safe_json
from mokume.core.logger import get_logger

logger = get_logger("mokume.reports.comparison")


def generate_comparison_report(
    workflows: list[dict],
    output_html: str = "workflow_comparison.html",
    title: str = "Workflow Comparison",
    marker_genes: Optional[dict[str, dict]] = None,
) -> str:
    """
    Generate interactive HTML comparing multiple quantification workflows.

    Parameters
    ----------
    workflows : list[dict]
        Each dict has keys:
        - name: str (e.g., "IRS_globalMedian")
        - protein_df: pd.DataFrame (wide format)
        - de_results: pd.DataFrame (optional)
        - sample_to_condition: dict
        - sample_to_plex: dict (optional)
        - is_log2: bool (default False)
    output_html : str
        Output path for the HTML report.
    title : str
        Report title.
    marker_genes : dict, optional
        Expected marker genes. Format: {accession: {"name": "GENE", "expected": "UP/DOWN"}}

    Returns
    -------
    str
        Path to the generated HTML file.
    """
    # Compute metrics for each workflow
    workflow_metrics = []
    for wf in workflows:
        metrics = compute_qc_metrics(
            protein_df=wf["protein_df"],
            sample_to_condition=wf["sample_to_condition"],
            sample_to_plex=wf.get("sample_to_plex"),
            de_results=wf.get("de_results"),
            is_log2=wf.get("is_log2", False),
        )
        workflow_metrics.append({
            "name": wf["name"],
            "metrics": metrics,
            "de_results": wf.get("de_results"),
        })

    # Compute cross-workflow concordance
    concordance = _compute_concordance(workflow_metrics)

    # Compute marker gene results
    marker_results = None
    if marker_genes:
        marker_results = _compute_marker_results(workflow_metrics, marker_genes)

    # Build HTML
    html = _build_comparison_html(
        title, workflow_metrics, concordance, marker_results, marker_genes,
    )

    with open(output_html, "w") as f:
        f.write(html)

    logger.info(f"Comparison report saved to {output_html}")
    return output_html


def _compute_concordance(workflow_metrics):
    """Compute pairwise DE concordance between workflows."""
    results = {}

    for i, wf_a in enumerate(workflow_metrics):
        for j, wf_b in enumerate(workflow_metrics):
            if j <= i:
                continue
            de_a = wf_a.get("de_results")
            de_b = wf_b.get("de_results")
            if de_a is None or de_b is None:
                continue

            key = f"{wf_a['name']} vs {wf_b['name']}"

            # DE overlap
            sig_a = set(de_a[de_a["significance"].isin(["UP", "DOWN"])]["ProteinName"])
            sig_b = set(de_b[de_b["significance"].isin(["UP", "DOWN"])]["ProteinName"])
            overlap = sig_a & sig_b
            only_a = sig_a - sig_b
            only_b = sig_b - sig_a

            # Log2FC correlation (all shared proteins)
            shared_proteins = set(de_a["ProteinName"]) & set(de_b["ProteinName"])
            if shared_proteins:
                merged = de_a[["ProteinName", "log2FC"]].merge(
                    de_b[["ProteinName", "log2FC"]],
                    on="ProteinName", suffixes=("_a", "_b"),
                )
                fc_corr = float(merged["log2FC_a"].corr(merged["log2FC_b"]))
                # Rank correlation
                fc_rank_corr = float(merged["log2FC_a"].corr(merged["log2FC_b"], method="spearman"))
                fc_data = {
                    "proteins": merged["ProteinName"].tolist(),
                    "fc_a": merged["log2FC_a"].tolist(),
                    "fc_b": merged["log2FC_b"].tolist(),
                }
            else:
                fc_corr = float("nan")
                fc_rank_corr = float("nan")
                fc_data = None

            # Concordance at top (CAT): for top-k DE genes by p-value,
            # what fraction overlap?
            cat_curve = _concordance_at_top(de_a, de_b)

            # Direction concordance: of shared DE genes, how many agree on direction?
            if overlap:
                merged_sig = de_a[de_a["ProteinName"].isin(overlap)][["ProteinName", "significance"]].merge(
                    de_b[de_b["ProteinName"].isin(overlap)][["ProteinName", "significance"]],
                    on="ProteinName", suffixes=("_a", "_b"),
                )
                direction_agree = (merged_sig["significance_a"] == merged_sig["significance_b"]).sum()
                direction_concordance = float(direction_agree / len(merged_sig))
            else:
                direction_concordance = float("nan")

            results[key] = {
                "n_a": len(sig_a),
                "n_b": len(sig_b),
                "overlap": len(overlap),
                "only_a": len(only_a),
                "only_b": len(only_b),
                "jaccard": float(len(overlap) / len(sig_a | sig_b)) if (sig_a | sig_b) else 0,
                "fc_pearson": fc_corr,
                "fc_spearman": fc_rank_corr,
                "direction_concordance": direction_concordance,
                "cat_curve": cat_curve,
                "fc_scatter": fc_data,
                "name_a": wf_a["name"],
                "name_b": wf_b["name"],
            }

    return results


def _concordance_at_top(de_a, de_b, max_k=200):
    """Concordance-At-Top: fraction of top-k DE genes shared."""
    # Sort by p-value
    sorted_a = de_a.sort_values("adj_pvalue")["ProteinName"].tolist()
    sorted_b = de_b.sort_values("adj_pvalue")["ProteinName"].tolist()

    max_k = min(max_k, len(sorted_a), len(sorted_b))
    ks = list(range(10, max_k + 1, 5))
    fractions = []

    for k in ks:
        top_a = set(sorted_a[:k])
        top_b = set(sorted_b[:k])
        frac = len(top_a & top_b) / k if k > 0 else 0
        fractions.append(frac)

    return {"k": ks, "fraction": fractions}


def _compute_marker_results(workflow_metrics, marker_genes):
    """Evaluate known marker genes across workflows."""
    results = {}
    for wf in workflow_metrics:
        de = wf.get("de_results")
        if de is None:
            continue
        wf_results = {}
        for acc, info in marker_genes.items():
            match = de[de["ProteinName"] == acc]
            if not match.empty:
                r = match.iloc[0]
                detected = r["significance"] in ["UP", "DOWN"]
                correct_direction = r["significance"] == info["expected"]
                wf_results[acc] = {
                    "gene": info["name"],
                    "expected": info["expected"],
                    "log2FC": float(r["log2FC"]),
                    "adj_pvalue": float(r["adj_pvalue"]),
                    "significance": str(r["significance"]),
                    "detected": detected,
                    "correct_direction": correct_direction,
                }
            else:
                wf_results[acc] = {
                    "gene": info["name"],
                    "expected": info["expected"],
                    "log2FC": None,
                    "adj_pvalue": None,
                    "significance": "NOT FOUND",
                    "detected": False,
                    "correct_direction": False,
                }
        results[wf["name"]] = wf_results
    return results


def _build_comparison_html(title, workflow_metrics, concordance, marker_results, marker_genes):
    """Build the comparison report HTML."""
    n_workflows = len(workflow_metrics)

    # Build summary table rows
    summary_rows = []
    for wf in workflow_metrics:
        m = wf["metrics"]
        s = m["summary"]
        sil = m["silhouette"]
        pvca = m["pvca"]
        de_q = m.get("de_quality", {})

        # Overall condition CV (median across conditions)
        all_cvs = [v["median_cv"] for v in m["condition_cv"].values()]
        median_cv = float(np.nanmedian(all_cvs)) if all_cvs else float("nan")

        # Overall intra-condition correlation
        all_corrs = list(m["intra_condition_corr"].values())
        median_corr = float(np.nanmedian(all_corrs)) if all_corrs else float("nan")

        summary_rows.append({
            "name": wf["name"],
            "n_proteins": s["n_proteins"],
            "silhouette": sil["score"],
            "var_condition": pvca["condition_pct"],
            "var_plex": pvca["plex_pct"],
            "median_cv": median_cv,
            "replicate_corr": median_corr,
            "n_de": de_q.get("n_de", 0),
            "n_up": de_q.get("n_up", 0),
            "n_down": de_q.get("n_down", 0),
            "pi1": de_q.get("pi1_estimate", float("nan")),
            "effect_size": de_q.get("effect_size", {}).get("median", 0),
            "direction_balance": de_q.get("direction_balance", float("nan")),
        })

    # Build PCA data per workflow
    pca_data = {}
    for wf in workflow_metrics:
        pca_data[wf["name"]] = wf["metrics"]["pca"]

    # Build concordance data
    concordance_json = {}
    for key, data in concordance.items():
        concordance_json[key] = {
            k: v for k, v in data.items()
            if k != "fc_scatter"  # exclude large data from main JSON
        }

    # FC scatter data (separate for size)
    fc_scatter_data = {}
    for key, data in concordance.items():
        if data.get("fc_scatter"):
            fc_scatter_data[key] = data["fc_scatter"]

    # Marker gene table
    marker_html = ""
    if marker_results and marker_genes:
        marker_html = _build_marker_table(marker_results, marker_genes, workflow_metrics)

    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{html.escape(title)}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
               background: #f5f5f5; color: #333; }}
        .header {{ background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
                   color: white; padding: 24px 30px; }}
        .header h1 {{ font-size: 24px; }}
        .header p {{ opacity: 0.8; margin-top: 6px; font-size: 13px; }}
        .container {{ max-width: 1500px; margin: 0 auto; padding: 20px; }}
        .card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                 padding: 16px; margin-bottom: 20px; }}
        .card h3 {{ font-size: 15px; color: #2c3e50; border-bottom: 2px solid #ecf0f1;
                    padding-bottom: 8px; margin-bottom: 12px; }}
        .card h4 {{ font-size: 13px; color: #7f8c8d; margin: 12px 0 8px; }}
        .row {{ display: flex; gap: 20px; margin-bottom: 20px; flex-wrap: wrap; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
        th {{ background: #34495e; color: white; padding: 8px 10px; text-align: center;
              font-size: 11px; }}
        td {{ padding: 6px 10px; border-bottom: 1px solid #eee; text-align: center; }}
        tr:hover {{ background: #f8f9fa; }}
        .best {{ background: #d5f5e3; font-weight: bold; }}
        .worst {{ background: #fadbd8; }}
        .good {{ color: #27ae60; }}
        .bad {{ color: #e74c3c; }}
        .neutral {{ color: #7f8c8d; }}
        .section {{ margin-bottom: 30px; }}
        .section-title {{ font-size: 18px; color: #2c3e50; margin-bottom: 16px;
                          padding-left: 10px; border-left: 4px solid #3498db; }}
        .pca-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
                     gap: 16px; }}
        .metric-explanation {{ font-size: 12px; color: #7f8c8d; margin-top: 8px;
                               line-height: 1.5; }}
    </style>
</head>
<body>

<div class="header">
    <h1>{html.escape(title)}</h1>
    <p>Comparing {n_workflows} quantification workflows. Higher silhouette, higher condition variance,
       lower plex variance, and lower CV indicate better performance.</p>
</div>

<div class="container">

    <!-- Section 1: Summary Table -->
    <div class="section">
        <h2 class="section-title">Summary Metrics</h2>
        <div class="card">
            <table id="summary-table">
                <thead><tr>
                    <th>Workflow</th>
                    <th>Proteins</th>
                    <th title="Silhouette score: how well samples cluster by condition (-1 to 1, higher=better)">Silhouette</th>
                    <th title="% variance explained by biological condition (higher=better)">% Var Condition</th>
                    <th title="% variance explained by technical plex (lower=better)">% Var Plex</th>
                    <th title="Median coefficient of variation within conditions (lower=better)">Median CV</th>
                    <th title="Median pairwise Pearson correlation within conditions (higher=better)">Replicate Corr</th>
                    <th>DE (UP/DOWN)</th>
                    <th title="Estimated proportion of true positives (Storey's pi1)">pi1</th>
                    <th title="Median |log2FC| of significant genes">Effect Size</th>
                    <th title="min(UP,DOWN)/max(UP,DOWN) - balance of direction (1=symmetric)">Balance</th>
                </tr></thead>
                <tbody id="summary-tbody"></tbody>
            </table>
            <p class="metric-explanation">
                <b>Silhouette</b>: Measures how well samples cluster by condition (1=perfect, 0=random, negative=wrong clusters).<br>
                <b>% Var Condition / Plex</b>: PVCA-style variance decomposition. Good workflow has high condition %, low plex %.<br>
                <b>Replicate Corr</b>: Median pairwise Pearson r between replicates within same condition (higher=more reproducible).<br>
                <b>pi1</b>: Estimated fraction of true positives from p-value distribution (Storey's method).<br>
                <b>Balance</b>: Ratio of smaller to larger DE direction (1=symmetric). Extreme asymmetry may indicate bias.
            </p>
        </div>
    </div>

    <!-- Section 2: PCA Comparison -->
    <div class="section">
        <h2 class="section-title">PCA by Condition + Plex</h2>
        <div class="pca-grid" id="pca-container"></div>
    </div>

    <!-- Section 3: DE Concordance -->
    <div class="section">
        <h2 class="section-title">DE Concordance</h2>
        <div class="row">
            <div style="flex:1">
                <div class="card">
                    <h3>Pairwise Overlap</h3>
                    <table id="concordance-table"></table>
                </div>
            </div>
            <div style="flex:1">
                <div class="card">
                    <h3>Log2FC Correlation (all shared proteins)</h3>
                    <div id="fc-scatter" style="height:400px"></div>
                </div>
            </div>
        </div>
        <div class="card">
            <h3>Concordance at Top</h3>
            <p class="metric-explanation" style="margin-bottom:12px">
                For the top-k most significant proteins (by p-value), what fraction overlap
                between methods? Plateau near 1.0 = strong agreement on most significant hits.
            </p>
            <div id="cat-plot" style="height:350px"></div>
        </div>
    </div>

    <!-- Section 4: Known Markers -->
    {f'''
    <div class="section">
        <h2 class="section-title">Known Marker Gene Panel</h2>
        <div class="card">
            <p class="metric-explanation" style="margin-bottom:12px">
                Known NASH-related genes. A good workflow should detect these with the expected
                direction and strong significance.
            </p>
            {marker_html}
        </div>
    </div>
    ''' if marker_html else ''}

</div>

<script>
const summaryRows = {_safe_json(summary_rows)};
const pcaData = {_safe_json(pca_data)};
const concordanceData = {_safe_json(concordance_json)};
const fcScatterData = {_safe_json(fc_scatter_data)};

const palette = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b'];
const plexSymbols = ['circle', 'square', 'diamond', 'cross', 'triangle-up', 'star'];

// ---- Summary Table ----
(function() {{
    const tbody = document.getElementById('summary-tbody');
    // Find best/worst for each metric
    const metrics = ['silhouette', 'var_condition', 'var_plex', 'median_cv', 'replicate_corr', 'pi1', 'effect_size'];
    const higher_is_better = ['silhouette', 'var_condition', 'replicate_corr', 'pi1', 'effect_size'];
    const lower_is_better = ['var_plex', 'median_cv'];

    const best = {{}};
    const worst = {{}};
    metrics.forEach(m => {{
        const vals = summaryRows.map(r => r[m]).filter(v => !isNaN(v) && v !== null);
        if (vals.length === 0) return;
        if (higher_is_better.includes(m)) {{
            best[m] = Math.max(...vals);
            worst[m] = Math.min(...vals);
        }} else if (lower_is_better.includes(m)) {{
            best[m] = Math.min(...vals);
            worst[m] = Math.max(...vals);
        }}
    }});

    function cellClass(metric, value) {{
        if (isNaN(value) || value === null) return '';
        if (value === best[metric]) return 'best';
        if (value === worst[metric] && summaryRows.length > 2) return 'worst';
        return '';
    }}

    summaryRows.forEach(row => {{
        const tr = document.createElement('tr');
        tr.innerHTML = `
            <td style="text-align:left;font-weight:bold">${{row.name}}</td>
            <td>${{row.n_proteins}}</td>
            <td class="${{cellClass('silhouette', row.silhouette)}}">${{row.silhouette.toFixed(3)}}</td>
            <td class="${{cellClass('var_condition', row.var_condition)}}">${{row.var_condition.toFixed(1)}}%</td>
            <td class="${{cellClass('var_plex', row.var_plex)}}">${{row.var_plex.toFixed(1)}}%</td>
            <td class="${{cellClass('median_cv', row.median_cv)}}">${{row.median_cv.toFixed(3)}}</td>
            <td class="${{cellClass('replicate_corr', row.replicate_corr)}}">${{row.replicate_corr.toFixed(3)}}</td>
            <td>${{row.n_de}} (${{row.n_up}}/${{row.n_down}})</td>
            <td class="${{cellClass('pi1', row.pi1)}}">${{isNaN(row.pi1) ? 'N/A' : row.pi1.toFixed(2)}}</td>
            <td class="${{cellClass('effect_size', row.effect_size)}}">${{row.effect_size.toFixed(2)}}</td>
            <td>${{isNaN(row.direction_balance) ? 'N/A' : row.direction_balance.toFixed(2)}}</td>
        `;
        tbody.appendChild(tr);
    }});
}})();

// ---- PCA Plots ----
(function() {{
    const container = document.getElementById('pca-container');
    Object.keys(pcaData).forEach((wfName, wfIdx) => {{
        const pca = pcaData[wfName];
        const div = document.createElement('div');
        div.className = 'card';
        div.innerHTML = `<div id="pca-${{wfIdx}}" style="height:380px"></div>`;
        container.appendChild(div);

        const conditions = [...new Set(pca.conditions)].sort();
        const plexes = [...new Set(pca.plexes)].sort();
        const condColors = {{}};
        const plexSymMap = {{}};
        conditions.forEach((c, i) => condColors[c] = palette[i % palette.length]);
        plexes.forEach((p, i) => plexSymMap[p] = plexSymbols[i % plexSymbols.length]);

        const traces = [];
        conditions.forEach(cond => {{
            plexes.forEach(plex => {{
                const idx = pca.samples.map((s, i) => i)
                    .filter(i => pca.conditions[i] === cond && pca.plexes[i] === plex);
                if (idx.length === 0) return;
                traces.push({{
                    type: 'scatter', mode: 'markers+text',
                    x: idx.map(i => pca.pc1[i]),
                    y: idx.map(i => pca.pc2[i]),
                    text: idx.map(i => pca.samples[i]),
                    textposition: 'top right',
                    textfont: {{ size: 7 }},
                    marker: {{
                        color: condColors[cond],
                        symbol: plexSymMap[plex],
                        size: 10,
                        line: {{ width: 1, color: '#333' }}
                    }},
                    name: cond + ' / ' + plex,
                }});
            }});
        }});

        Plotly.newPlot(`pca-${{wfIdx}}`, traces, {{
            title: {{ text: wfName, font: {{ size: 14 }} }},
            xaxis: {{ title: 'PC1 (' + pca.var_explained[0].toFixed(1) + '%)' }},
            yaxis: {{ title: 'PC2 (' + pca.var_explained[1].toFixed(1) + '%)' }},
            legend: {{ font: {{ size: 9 }}, orientation: 'h', y: -0.15 }},
            margin: {{ t: 50, b: 80, l: 60, r: 20 }},
        }}, {{responsive: true}});
    }});
}})();

// ---- Concordance Table ----
(function() {{
    const table = document.getElementById('concordance-table');
    let html = '<thead><tr><th>Comparison</th><th>DE A</th><th>DE B</th><th>Overlap</th>' +
               '<th>Jaccard</th><th>FC Pearson</th><th>FC Spearman</th><th>Direction</th></tr></thead><tbody>';

    Object.keys(concordanceData).forEach(key => {{
        const d = concordanceData[key];
        html += `<tr>
            <td style="text-align:left">${{key}}</td>
            <td>${{d.n_a}}</td><td>${{d.n_b}}</td>
            <td><b>${{d.overlap}}</b></td>
            <td>${{d.jaccard.toFixed(3)}}</td>
            <td>${{d.fc_pearson.toFixed(3)}}</td>
            <td>${{d.fc_spearman.toFixed(3)}}</td>
            <td>${{isNaN(d.direction_concordance) ? 'N/A' : (d.direction_concordance * 100).toFixed(1) + '%'}}</td>
        </tr>`;
    }});
    html += '</tbody>';
    table.innerHTML = html;
}})();

// ---- FC Scatter ----
(function() {{
    const keys = Object.keys(fcScatterData);
    if (keys.length === 0) return;
    // Show first pair
    const key = keys[0];
    const d = fcScatterData[key];
    const cd = concordanceData[key];

    Plotly.newPlot('fc-scatter', [{{
        type: 'scatter', mode: 'markers',
        x: d.fc_a, y: d.fc_b,
        text: d.proteins,
        marker: {{ size: 4, color: '#3498db', opacity: 0.5 }},
        hovertemplate: '<b>%{{text}}</b><br>%{{xaxis.title.text}}: %{{x:.2f}}<br>%{{yaxis.title.text}}: %{{y:.2f}}<extra></extra>',
    }}, {{
        type: 'scatter', mode: 'lines',
        x: [-10, 10], y: [-10, 10],
        line: {{ dash: 'dash', color: 'red', width: 1 }},
        showlegend: false,
    }}], {{
        title: {{ text: `log2FC: ${{cd.name_a}} vs ${{cd.name_b}} (r=${{cd.fc_pearson.toFixed(3)}})`, font: {{ size: 13 }} }},
        xaxis: {{ title: cd.name_a + ' log2FC', range: [-8, 8] }},
        yaxis: {{ title: cd.name_b + ' log2FC', range: [-8, 8] }},
        margin: {{ t: 50, b: 50, l: 60, r: 20 }},
    }}, {{responsive: true}});
}})();

// ---- CAT Curve ----
(function() {{
    const traces = [];
    Object.keys(concordanceData).forEach((key, i) => {{
        const cat = concordanceData[key].cat_curve;
        if (!cat) return;
        traces.push({{
            type: 'scatter', mode: 'lines',
            x: cat.k, y: cat.fraction,
            name: key,
            line: {{ width: 2, color: palette[i % palette.length] }},
        }});
    }});
    if (traces.length === 0) return;
    Plotly.newPlot('cat-plot', traces, {{
        title: {{ text: 'Concordance at Top-k', font: {{ size: 14 }} }},
        xaxis: {{ title: 'Top k proteins (by p-value)' }},
        yaxis: {{ title: 'Fraction overlap', range: [0, 1] }},
        margin: {{ t: 50, b: 50, l: 60, r: 20 }},
    }}, {{responsive: true}});
}})();
</script>

<div style="text-align:center; padding:20px; color:#999; font-size:11px;">
    Generated by <a href="https://github.com/bigbio/mokume" style="color:#666;">mokume</a>
</div>
</body>
</html>"""


def _build_marker_table(marker_results, marker_genes, workflow_metrics):
    """Build HTML table for marker gene results."""
    wf_names = [wf["name"] for wf in workflow_metrics if wf["name"] in marker_results]
    if not wf_names:
        return ""

    # Header row
    header = "<tr><th>Gene</th><th>Accession</th><th>Expected</th>"
    for wf in wf_names:
        header += f"<th colspan='2'>{wf}</th>"
    header += "</tr>"

    # Sub-header
    subheader = "<tr><th></th><th></th><th></th>"
    for _ in wf_names:
        subheader += "<th>log2FC</th><th>Significance</th>"
    subheader += "</tr>"

    rows = ""
    for acc, info in marker_genes.items():
        row = f"<tr><td><b>{info['name']}</b></td><td>{acc}</td>"
        row += f"<td class='{'good' if info['expected'] == 'UP' else 'bad'}'>{info['expected']}</td>"

        for wf in wf_names:
            r = marker_results[wf].get(acc, {})
            if r.get("log2FC") is not None:
                fc = f"{r['log2FC']:.2f}"
                sig = r["significance"]
                css = "good" if r.get("correct_direction") else ("bad" if r.get("detected") else "neutral")
                row += f"<td>{fc}</td><td class='{css}'>{sig}</td>"
            else:
                row += "<td>-</td><td class='neutral'>NOT FOUND</td>"

        row += "</tr>"
        rows += row

    return f"""
    <table>
        <thead>{header}{subheader}</thead>
        <tbody>{rows}</tbody>
    </table>
    """
