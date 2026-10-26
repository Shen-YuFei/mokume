"""Shared palettes and presentation constants for static plots and reports."""

CONDITION_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#56B4E9",
    "#CC79A7",
)
CONDITION_MARKERS = ("o", "s", "^", "D", "v", "P")
CONDITION_SYMBOLS = (
    "circle",
    "square",
    "triangle-up",
    "diamond",
    "triangle-down",
    "cross",
)
# Past six groups, cycle colors before shapes so neighbouring groups never share a hue.
CONDITION_STYLE_INDICES = tuple((index, index) for index in range(6)) + tuple(
    (color, (color + shift) % 6) for shift in range(1, 6) for color in range(6)
)
PLEX_SYMBOLS = tuple(range(53))
SIGNIFICANCE_COLORS = {
    "UP": "#D55E00",
    "DOWN": "#0072B2",
    "Unchanged": "#BBBBBB",
    "NotTested": "#777777",
}
MISSING_COLOR = "#BDBDBD"
# Single-hue blue ramp, light (low) to dark (high), for magnitude heatmaps.
SEQUENTIAL_COLORS = (
    "#cde2fb",
    "#9ec5f4",
    "#6da7ec",
    "#3987e5",
    "#256abf",
    "#184f95",
    "#0d366b",
)
SINGLE_SERIES_COLOR = "#0072B2"
SAMPLE_POINT_COLOR = "#52514e"
TEXT_COLOR = "#333333"
GRID_COLOR = "#E5E5E5"
AXIS_COLOR = "#D9D9D9"
WHITE = "#ffffff"
PEPTIDE_COLOR = "#D0D0D0"
INACTIVE_COLOR = "#E5E5E5"
PLOTLY_AXIS_STYLE = {
    "showgrid": True,
    "gridcolor": GRID_COLOR,
    "gridwidth": 1,
    "griddash": "solid",
    "showline": True,
    "linecolor": AXIS_COLOR,
    "linewidth": 1,
    "zeroline": False,
}
PLOTLY_TEMPLATE = {
    "layout": {
        "font": {"color": TEXT_COLOR},
        "hoverlabel": {
            "bgcolor": "#ffffff",
            "bordercolor": AXIS_COLOR,
            "font": {"color": TEXT_COLOR},
        },
        "xaxis": PLOTLY_AXIS_STYLE,
        "yaxis": PLOTLY_AXIS_STYLE,
        "bargap": 0.3,
        "paper_bgcolor": "#ffffff",
        "plot_bgcolor": "#ffffff",
    },
    "data": {
        "scatter": [
            {
                "marker": {"size": 8, "line": {"width": 2, "color": "#ffffff"}},
                "textfont": {"color": TEXT_COLOR},
            }
        ],
        "bar": [{"marker": {"line": {"width": 0}}, "textfont": {"color": TEXT_COLOR}}],
    },
}

REPORT_CSS = f"""
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       background: #f5f5f5; color: #333; }}
.header {{ background: #2c3e50; color: #ffffff; padding: 20px 30px; }}
.header h1 {{ font-size: 24px; margin-bottom: 8px; }}
.stat-box, .header-stat {{ background: #f5f5f5; color: {TEXT_COLOR}; }}
.card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
         padding: 16px; margin-bottom: 20px; min-width: 0; }}
.card h3 {{ margin-bottom: 12px; font-size: 15px; color: {TEXT_COLOR};
           border-bottom: 2px solid #ecf0f1; padding-bottom: 8px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th {{ background: #ecf0f1; color: {TEXT_COLOR}; padding: 8px 10px; text-align: left; }}
td {{ padding: 6px 10px; border-bottom: 1px solid #eee; }}
tr:hover {{ background: #f8f9fa; }}
.sig-up, .sig-down {{ font-weight: bold; }}
.sig-up::before, .sig-down::before, .sig-unchanged::before, .sig-nottested::before {{
    content: ''; display: inline-block; width: 0.7em; height: 0.7em; margin-right: 0.4em;
}}
.sig-up::before {{ background: {SIGNIFICANCE_COLORS["UP"]}; }}
.sig-down::before {{ background: {SIGNIFICANCE_COLORS["DOWN"]}; }}
.sig-unchanged::before {{ background: {SIGNIFICANCE_COLORS["Unchanged"]}; }}
.sig-nottested::before {{ background: {SIGNIFICANCE_COLORS["NotTested"]}; }}
.stat-up {{ border-left: 3px solid {SIGNIFICANCE_COLORS["UP"]}; }}
.stat-down {{ border-left: 3px solid {SIGNIFICANCE_COLORS["DOWN"]}; }}
.info-text {{ font-size: 12px; color: #666; margin-top: 8px; }}
.plot-grid {{ display: grid; gap: 0 20px;
              grid-template-columns: repeat(auto-fit, minmax(min(100%, 560px), 1fr)); }}
.plot-grid .wide {{ grid-column: 1 / -1; }}
"""
