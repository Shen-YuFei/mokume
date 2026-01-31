#!/usr/bin/env python3
"""
Generate a workflow diagram for the mokume protein quantification library.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


def create_workflow_diagram():
    """Create a clean workflow diagram for mokume."""

    fig, ax = plt.subplots(1, 1, figsize=(16, 11))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 11)
    ax.axis('off')
    ax.set_facecolor('white')

    # Color palette
    colors = {
        'input': '#FFECB3',
        'input_border': '#FFA000',
        'process': '#BBDEFB',
        'process_border': '#1976D2',
        'norm': '#F8BBD0',
        'norm_border': '#C2185B',
        'norm_item': '#FCE4EC',
        'quant': '#C8E6C9',
        'quant_border': '#388E3C',
        'quant_item': '#E8F5E9',
        'output': '#D1C4E9',
        'output_border': '#7B1FA2',
        'dlfq': '#E1F5FE',
        'dlfq_border': '#0288D1',
        'dlfq_item': '#B3E5FC',
        'arrow': '#455A64',
    }

    def draw_box(x, y, w, h, fc, ec, lw=2):
        """Draw a rounded box."""
        box = FancyBboxPatch((x, y), w, h,
                            boxstyle="round,pad=0.02,rounding_size=0.1",
                            facecolor=fc, edgecolor=ec, linewidth=lw)
        ax.add_patch(box)

    def draw_arrow(x1, y1, x2, y2, style='-', color=None):
        """Draw arrow."""
        if color is None:
            color = colors['arrow']
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='-|>', color=color, lw=1.5,
                                 linestyle=style))

    # ==================== TITLE ====================
    ax.text(8, 10.5, 'mokume', ha='center', va='center',
           fontsize=26, fontweight='bold', color='#1565C0')
    ax.text(8, 10.0, 'Protein Quantification Workflow', ha='center', va='center',
           fontsize=13, color='#666666')

    # ==================== ROW 1: INPUT FILES ====================
    y_row1 = 7.0

    # Feature file
    draw_box(0.5, y_row1, 1.4, 1.1, colors['input'], colors['input_border'])
    ax.text(1.2, y_row1 + 0.7, 'Feature', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['input_border'])
    ax.text(1.2, y_row1 + 0.35, 'Parquet', ha='center', va='center',
           fontsize=8, color='#666666')

    # SDRF file
    draw_box(0.5, y_row1 - 1.8, 1.4, 1.1, colors['input'], colors['input_border'])
    ax.text(1.2, y_row1 - 1.1, 'SDRF', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['input_border'])
    ax.text(1.2, y_row1 - 1.45, 'TSV', ha='center', va='center',
           fontsize=8, color='#666666')

    # ==================== STAGE 1: features2peptides ====================
    draw_box(2.5, y_row1 - 0.8, 1.8, 1.3, colors['process'], colors['process_border'])
    ax.text(3.4, y_row1 - 0.05, 'features2', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])
    ax.text(3.4, y_row1 - 0.4, 'peptides', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])

    draw_arrow(1.9, y_row1 + 0.5, 2.5, y_row1 + 0.1)
    draw_arrow(1.9, y_row1 - 1.2, 2.5, y_row1 - 0.4)

    # ==================== FEATURE NORMALIZATION ====================
    norm1_x, norm1_y = 5.0, 7.5
    norm1_w, norm1_h = 2.0, 2.2

    draw_box(norm1_x, norm1_y, norm1_w, norm1_h, colors['norm'], colors['norm_border'], lw=1.5)
    ax.text(norm1_x + norm1_w/2, norm1_y + norm1_h - 0.3, 'Feature Norm',
           ha='center', va='center', fontsize=9, fontweight='bold',
           color=colors['norm_border'], fontstyle='italic')

    norm1_items = ['Median', 'Mean', 'IQR', 'Global']
    item_h = 0.38
    for i, item in enumerate(norm1_items):
        iy = norm1_y + norm1_h - 0.75 - i * (item_h + 0.05)
        draw_box(norm1_x + 0.15, iy, norm1_w - 0.3, item_h, colors['norm_item'], colors['norm_border'], lw=1)
        ax.text(norm1_x + norm1_w/2, iy + item_h/2, item, ha='center', va='center',
               fontsize=8, color='#333333')

    # ==================== PEPTIDE NORMALIZATION ====================
    norm2_x, norm2_y = 5.0, 4.6
    norm2_w, norm2_h = 2.0, 1.9

    draw_box(norm2_x, norm2_y, norm2_w, norm2_h, colors['norm'], colors['norm_border'], lw=1.5)
    ax.text(norm2_x + norm2_w/2, norm2_y + norm2_h - 0.3, 'Peptide Norm',
           ha='center', va='center', fontsize=9, fontweight='bold',
           color=colors['norm_border'], fontstyle='italic')

    norm2_items = ['GlobalMedian', 'ConditionMedian', 'Log2']
    for i, item in enumerate(norm2_items):
        iy = norm2_y + norm2_h - 0.75 - i * (item_h + 0.05)
        draw_box(norm2_x + 0.15, iy, norm2_w - 0.3, item_h, colors['norm_item'], colors['norm_border'], lw=1)
        ax.text(norm2_x + norm2_w/2, iy + item_h/2, item, ha='center', va='center',
               fontsize=8, color='#333333')

    draw_arrow(4.3, y_row1 + 0.0, 5.0, 8.2)
    draw_arrow(4.3, y_row1 - 0.5, 5.0, 5.7)

    # ==================== STAGE 2: Peptide Summary ====================
    draw_box(7.5, y_row1 - 0.8, 1.8, 1.3, colors['process'], colors['process_border'])
    ax.text(8.4, y_row1 - 0.05, 'Peptide', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])
    ax.text(8.4, y_row1 - 0.4, 'Summary', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])

    draw_arrow(7.0, 8.2, 7.5, y_row1 + 0.1)
    draw_arrow(7.0, 5.7, 7.5, y_row1 - 0.4)

    # ==================== STAGE 3: peptides2protein ====================
    draw_box(10.0, y_row1 - 0.8, 1.8, 1.3, colors['process'], colors['process_border'])
    ax.text(10.9, y_row1 - 0.05, 'peptides2', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])
    ax.text(10.9, y_row1 - 0.4, 'protein', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['process_border'])

    draw_arrow(9.3, y_row1 - 0.2, 10.0, y_row1 - 0.2)

    # ==================== QUANTIFICATION METHODS ====================
    quant_x, quant_y = 12.5, 4.2
    quant_w, quant_h = 3.0, 5.5

    draw_box(quant_x, quant_y, quant_w, quant_h, colors['quant'], colors['quant_border'], lw=2)
    ax.text(quant_x + quant_w/2, quant_y + quant_h - 0.4, 'Quantification',
           ha='center', va='center', fontsize=12, fontweight='bold',
           color=colors['quant_border'])

    methods = ['iBAQ', 'DirectLFQ', 'MaxLFQ', 'Top3', 'TopN', 'Sum']
    method_h = 0.6
    method_start = quant_y + quant_h - 1.1

    for i, method in enumerate(methods):
        my = method_start - i * (method_h + 0.15)
        draw_box(quant_x + 0.25, my, quant_w - 0.5, method_h, colors['quant_item'], colors['quant_border'], lw=1.2)
        ax.text(quant_x + quant_w/2, my + method_h/2, method,
               ha='center', va='center', fontsize=10, fontweight='bold',
               color=colors['quant_border'])

    draw_arrow(11.8, y_row1 - 0.2, 12.5, 8.0)

    # ==================== DIRECTLFQ NORMALIZATION ====================
    dlfq_x, dlfq_y = 9.8, 2.5
    dlfq_w, dlfq_h = 2.2, 2.2

    draw_box(dlfq_x, dlfq_y, dlfq_w, dlfq_h, colors['dlfq'], colors['dlfq_border'], lw=1.5)
    ax.text(dlfq_x + dlfq_w/2, dlfq_y + dlfq_h - 0.3, 'DirectLFQ Norm',
           ha='center', va='center', fontsize=9, fontweight='bold',
           color=colors['dlfq_border'], fontstyle='italic')

    dlfq_items = ['Hierarchical', 'Variance-guided', 'Pairwise align.']
    dlfq_item_h = 0.45
    for i, item in enumerate(dlfq_items):
        iy = dlfq_y + dlfq_h - 0.75 - i * (dlfq_item_h + 0.08)
        draw_box(dlfq_x + 0.15, iy, dlfq_w - 0.3, dlfq_item_h, colors['dlfq_item'], colors['dlfq_border'], lw=1)
        ax.text(dlfq_x + dlfq_w/2, iy + dlfq_item_h/2, item, ha='center', va='center',
               fontsize=8, color='#333333')

    draw_arrow(12.0, 3.6, 12.5, 5.8, style='--', color=colors['dlfq_border'])

    # ==================== FASTA & ENZYMATIC DIGESTION ====================
    # FASTA file
    draw_box(0.5, 1.2, 1.4, 1.1, colors['input'], colors['input_border'])
    ax.text(1.2, 1.85, 'FASTA', ha='center', va='center',
           fontsize=10, fontweight='bold', color=colors['input_border'])
    ax.text(1.2, 1.5, 'Protein DB', ha='center', va='center',
           fontsize=8, color='#666666')

    draw_box(3.2, 1.2, 2.2, 1.1, colors['process'], colors['process_border'])
    ax.text(4.3, 1.85, 'Enzymatic', ha='center', va='center',
           fontsize=9, fontweight='bold', color=colors['process_border'])
    ax.text(4.3, 1.5, 'Digestion', ha='center', va='center',
           fontsize=9, fontweight='bold', color=colors['process_border'])

    draw_arrow(1.9, 1.75, 3.2, 1.75)
    draw_arrow(5.4, 1.75, 12.5, 4.8, style='--', color='#1976D2')

    ax.text(8.8, 2.0, 'Required for iBAQ', ha='center', va='center',
           fontsize=8, fontstyle='italic', color='#1976D2')

    # ==================== OUTPUT ====================
    draw_box(13.0, 0.8, 2.2, 1.5, colors['output'], colors['output_border'])
    ax.text(14.1, 1.8, 'Report Matrix', ha='center', va='center',
           fontsize=11, fontweight='bold', color=colors['output_border'])
    ax.text(14.1, 1.35, 'TSV / Parquet', ha='center', va='center',
           fontsize=8, color='#666666')
    ax.text(14.1, 1.05, 'AnnData', ha='center', va='center',
           fontsize=8, color='#666666')

    draw_arrow(14.0, 4.2, 14.0, 2.3)

    # ==================== LEGEND ====================
    legend_y = 0.2
    legend_items = [
        (0.5, colors['input'], colors['input_border'], 'Input'),
        (2.2, colors['process'], colors['process_border'], 'Process'),
        (4.2, colors['norm'], colors['norm_border'], 'Normalization'),
        (6.5, colors['quant'], colors['quant_border'], 'Quantification'),
        (8.8, colors['output'], colors['output_border'], 'Output'),
    ]

    for lx, fc, ec, label in legend_items:
        draw_box(lx, legend_y, 0.3, 0.3, fc, ec, lw=1.5)
        ax.text(lx + 0.45, legend_y + 0.15, label, ha='left', va='center',
               fontsize=9, color='#333333')

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    fig = create_workflow_diagram()

    fig.savefig('plots/mokume_workflow.png', dpi=250, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    fig.savefig('plots/mokume_workflow.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')

    print("Workflow diagram saved to plots/mokume_workflow.png and plots/mokume_workflow.pdf")
    plt.close()
