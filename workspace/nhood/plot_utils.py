


from glob import glob
import os
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from utils import load_patients_df
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.gridspec as gridspec


def plot_clustered_heatmap(all_stats, treatment, outfile, label_colors=None, legend_labels=None, group_prefixes=None):
    """
    Create a clustered heatmap with dendrogram on right, clustering combinations (y-axis).

    Args:
        label_colors: dict with group prefix keys specifying colors for each group.
                     Default is {'group1': 'green', 'group2': 'red'}
        legend_labels: dict with group prefix keys specifying legend labels for each group.
                      Default is {'group1': 'Responder', 'group2': 'Non-Responder'}
        group_prefixes: dict with 'group1' and 'group2' keys specifying the prefixes to use.
                       Default is {'group1': 'R_', 'group2': 'NR_'}
    """
    if group_prefixes is None:
        group_prefixes = {'group1': 'R_', 'group2': 'NR_'}
    if label_colors is None:
        label_colors = {'group1': 'green', 'group2': 'red'}
    if legend_labels is None:
        legend_labels = {'group1': 'Responder', 'group2': 'Non-Responder'}
    # Create figure with 4 regions
    fig = plt.figure(figsize=(4, 13))

    regions = ["Tumor", "Stroma"]

    for idx, region in enumerate(regions):
        sub_stats = all_stats[all_stats['Region'] == region]
        region_df = []
        for _, row in sub_stats.iterrows():
            responder_dict = row['Responder dict']
            non_responder_dict = row['Non-responder dict']
            region_df.append({
                'Combination': row['Combination'],
                **{f"{group_prefixes['group1']}{k}": v for k, v in responder_dict.items()},
                **{f"{group_prefixes['group2']}{k}": v for k, v in non_responder_dict.items()},
                'P-value': row['P-value']
            })
        region_df = pd.DataFrame(region_df)

        # Skip if no data
        if region_df.empty:
            continue

        # Extract p-values and combination names
        p_values = region_df['P-value'].values
        combinations = region_df['Combination'].values

        # Prepare heatmap data: set Combination as index, exclude P-value column
        heatmap_data = region_df.set_index('Combination').drop(columns=['P-value'])

        # Skip if no data
        if heatmap_data.empty:
            continue

        # Apply z-score normalization row-wise, preserving NaN for missing data
        # nan_policy='omit' calculates zscore only on non-NaN values
        def zscore_with_nan(row):
            """Apply zscore while preserving NaN values"""
            mask = ~np.isnan(row)
            if mask.sum() < 2:  # Need at least 2 values for zscore
                return row * 0  # Return zeros if insufficient data
            result = row.copy()
            valid_values = row[mask]
            if valid_values.std() == 0:  # No variation
                result[mask] = 0
            else:
                result[mask] = (valid_values - valid_values.mean()) / valid_values.std()
            return result

        def val_with_nan(row):
            """Same as the above function but without normalization"""
            mask = ~np.isnan(row)
            if mask.sum() < 2:  # Need at least 2 values for zscore
                return row * 0  # Return zeros if insufficient data
            result = row.copy()
            valid_values = row[mask]
            result[mask] = valid_values
            return result

        heatmap_data_zscore = heatmap_data.apply(zscore_with_nan, axis=1)

        # For clustering, fill NaN with 0 (neutral zscore)
        heatmap_data_zscore_for_clustering = heatmap_data_zscore.fillna(0)

        # Perform hierarchical clustering on combinations (rows)
        linkage_matrix = linkage(heatmap_data_zscore_for_clustering, method='average', metric='euclidean')

        # Create a gridspec for this region: heatmap, colorbar (no dendrogram)
        # Adjust spacing to prevent overlap
        height_per_region = 0.15
        gap_between_regions = 0.15
        top_position = 0.97 - idx * (height_per_region + gap_between_regions)
        bottom_position = top_position - height_per_region

        gs = gridspec.GridSpec(1, 2, figure=fig,
                               width_ratios=[1, 0.03],
                               wspace=0.08,  # Space for p-value labels
                               left=0.05,
                               right=0.95,
                               top=top_position,
                               bottom=bottom_position
                               )

        # Create heatmap axis (on the left)
        ax_heatmap = fig.add_subplot(gs[0])

        # Create dendrogram axis (in the middle)
        # ax_dendro = fig.add_subplot(gs[1])

        # Create colorbar axis (on the right)
        # ax_cbar = fig.add_subplot(gs[1])

        # Reorder rows based on dendrogram first (before plotting)
        dendro = dendrogram(linkage_matrix, no_plot=True)
        row_order = dendro['leaves']
        heatmap_data_clustered = heatmap_data_zscore.iloc[row_order, :]
        p_values_clustered = p_values[row_order]

        # Create the heatmap with z-score normalized and clustered data (no colorbar)
        im = ax_heatmap.pcolormesh(heatmap_data_clustered, cmap='bwr', vmin=-3, vmax=3)
        ax_heatmap.set_ylim(0, len(heatmap_data_clustered))
        ax_heatmap.set_xlim(0, len(heatmap_data_clustered.columns))
        ax_heatmap.invert_yaxis()

        # Set ticks for heatmap
        ax_heatmap.set_xticks(np.arange(len(heatmap_data_clustered.columns)) + 0.5)
        ax_heatmap.set_yticks(np.arange(len(heatmap_data_clustered)) + 0.5)

        # Set labels
        ax_heatmap.set_xticklabels(heatmap_data_clustered.columns, rotation=45, ha='right', fontsize=8)
        yticklabels = [label.replace('_', ' ') for label in heatmap_data_clustered.index]
        ax_heatmap.set_yticklabels(yticklabels, fontsize=8)

        # Now plot the dendrogram with the same clustering
        # dendro_plot = dendrogram(linkage_matrix, ax=ax_dendro, orientation='right',
        #                         no_labels=True, color_threshold=0, above_threshold_color='black')

        # Invert the dendrogram y-axis to match the heatmap orientation
        # ax_dendro.invert_yaxis()

        # ax_dendro.set_xticks([])
        # ax_dendro.set_yticks([])
        # ax_dendro.spines['top'].set_visible(False)
        # ax_dendro.spines['right'].set_visible(False)
        # ax_dendro.spines['bottom'].set_visible(False)
        # ax_dendro.spines['left'].set_visible(False)

        # Color x-axis labels based on group prefix
        labels = ax_heatmap.get_xticklabels()
        for label in labels:
            sample_name = label.get_text()
            if sample_name.startswith(group_prefixes['group1']):
                label.set_color(label_colors['group1'])
            elif sample_name.startswith(group_prefixes['group2']):
                label.set_color(label_colors['group2'])

        # Add colorbar at the rightmost position
        #cbar = plt.colorbar(im, cax=ax_cbar, )
        #cbar.set_label('Z-score', fontsize=9)

        # Add p-values on the right side for each combination (in clustered order)
        for i, p_val in enumerate(p_values_clustered):
            if pd.notna(p_val):
                color = 'red' if p_val < 0.05 else 'black'
                x_pos = len(heatmap_data_clustered.columns) + 0.1
                ax_heatmap.text(x_pos, i + 0.5, f"{p_val:.2f}", ha='left', va='center',
                       fontsize=8, transform=ax_heatmap.transData, color=color, clip_on=False)

        # Set title and labels
        # ax_dendro.set_title(f'{region} - {treatment} (Clustered)', fontsize=10, fontweight='bold',
        #                    loc='left', pad=10)
        ax_heatmap.set_title(f'{region} - {treatment}', fontsize=10, fontweight='bold',
                           loc='left', pad=10)
        ax_heatmap.set_ylabel('Cell Type Combinations', fontsize=9)
        ax_heatmap.set_xlabel('Samples', fontsize=9)

    # Create legend for the entire figure
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor=label_colors['group1'],
               markersize=8, label=legend_labels['group1']),
        Line2D([0], [0], marker='s', color='w', markerfacecolor=label_colors['group2'],
               markersize=8, label=legend_labels['group2'])
    ]
    #fig.legend(handles=legend_elements, loc='lower left', fontsize=10)

    # Save the figure
    plt.savefig(outfile.replace(".png", "_.png"), bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved clustered heatmap: {outfile}")


# THe same as in the script 601
def create_combined_plot_v2(
        per_sample_abundance,
        motif_association_df,
        per_motif_composition,
        response_data,
        outfile,
        treatment=None,
        m=None,
        mode='raw',
        prepare_motif_labels=True,
        fig_width=12
):
    # Create figure with gridspec to accommodate dendrogram
    from matplotlib.gridspec import GridSpec
    # compute height based on number of motifs
    n_motifs = per_motif_composition.shape[0]
    fig_height = max(4, n_motifs * 0.2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[4, 3], height_ratios=[1, 7, 0.5],
                  hspace=0.02, wspace=0.6)

    # Create axes
    ax1 = fig.add_subplot(gs[1, 0])  # Left heatmap (spans middle row)
    ax_dendro = fig.add_subplot(gs[0, 1])  # Dendrogram (top right)
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax_dendro)  # Right heatmap (middle right)

    # Left figure: abundance per sample
    # Calculate max absolute value for symmetric colorbar
    max_abs_val = np.abs(per_sample_abundance.values).max()
    sns.heatmap(
        per_sample_abundance.T, annot=False, cmap='bwr', center=0,
        yticklabels=1,
        vmin=-max_abs_val, vmax=max_abs_val, ax=ax1,
        cbar_kws={'label': 'Z-score', 'shrink': 0.3, 'pad': 0.13})
    ax1.set_title('Motif Abundance per Sample (Z-score normalized)')
    # Force all ticks to appear
    ax1.set_xticks(np.arange(per_sample_abundance.T.shape[1]) + 0.5)
    ax1.set_xticklabels(per_sample_abundance.T.columns, rotation=90)

    # Add p-values to the right of the heatmap
    if motif_association_df is not None:
        for i in range(len(motif_association_df)):
            y_pos = i
            x_pos = len(per_sample_abundance) + 0.3
            p_val = motif_association_df.iloc[i]['p_value']
            ax1.text(x_pos, y_pos + 0.5, f"{p_val:.3f}", ha='left', va='center', fontsize=10,
                    transform=ax1.transData, color='red' if p_val < 0.05 else 'black')

    # Rename motifs to CN1, CN2, etc. based on actual motif IDs
    if prepare_motif_labels:
        motif_labels = [f"CN{motif+1}" for motif in per_sample_abundance.columns]
    else:
        motif_labels = per_sample_abundance.columns.tolist()
    ax1.set_yticklabels(motif_labels, rotation=0)

    # Rotate x-ticks and set alignment
    ax1.tick_params(axis='x', rotation=45)
    plt.setp(ax1.get_xticklabels(), ha='right')

    # Color sample labels by response
    if response_data is not None:
        labels = ax1.get_xticklabels()
        for label in labels:
            sample_name = label.get_text().split("_(")[0]
            response = response_data.loc[sample_name]
            if response == 'pos':
                label.set_color('green')
            elif response == 'neg':
                label.set_color('red')

    # Right figure: motif composition
    colors = ['white', '#FF0000']
    cmap = LinearSegmentedColormap.from_list('white_to_red', colors)
    cbar_label = 'Log count' if mode == 'log' else 'Raw'
    if mode == 'log':
        per_motif_composition = np.log2(1 + per_motif_composition)
    sns.heatmap(per_motif_composition, annot=False, vmin=0, cmap=cmap,
                xticklabels=1, yticklabels=1,
                cbar_kws={'label': cbar_label, 'shrink': 0.3, 'pad': 0.028}, ax=ax2)
    ax2.set_title(f'Motif Phenotype Composition',  pad=6)

    # Rename motifs to CN1, CN2, etc. on right heatmap
    ax2.set_yticklabels(motif_labels, rotation=0)

    # Rotate x-ticks on right heatmap and replace underscores with spaces
    ax2.tick_params(axis='x', rotation=45)
    xticklabels = ax2.get_xticklabels()
    ax2.set_xticklabels([label.get_text().replace('_', ' ') for label in xticklabels], rotation=45, ha='right')

    # Hide the empty dendrogram axes
    ax_dendro.axis('off')

    # Add overall figure title
    if treatment and m:
        fig.suptitle(f'Cellular Neighborhood Analysis - {treatment} Treatment (m={m})',
                    fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, top=0.96)
    plt.savefig(outfile, bbox_inches='tight', dpi=300)


def create_combined_plot_v2_without_group(
        per_sample_abundance,
        motif_association_df,
        per_motif_composition,
        outfile,
        treatment=None,
        m=None,
        mode='raw',
        prepare_motif_labels=True,
        fig_width=12
):
    # Create figure with gridspec to accommodate dendrogram
    from matplotlib.gridspec import GridSpec
    # compute height based on number of motifs
    n_motifs = per_motif_composition.shape[0]
    fig_height = max(4, n_motifs * 0.2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[4, 3], height_ratios=[1, 7, 0.5],
                  hspace=0.02, wspace=0.6)

    # Create axes
    ax1 = fig.add_subplot(gs[1, 0])  # Left heatmap (spans middle row)
    ax_dendro = fig.add_subplot(gs[0, 1])  # Dendrogram (top right)
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax_dendro)  # Right heatmap (middle right)

    # Left figure: abundance per sample
    # Calculate max absolute value for symmetric colorbar
    max_abs_val = np.abs(per_sample_abundance.values).max()
    sns.heatmap(
        per_sample_abundance.T, annot=False, cmap='bwr', center=0,
        yticklabels=1,
        vmin=-max_abs_val, vmax=max_abs_val, ax=ax1,
        cbar_kws={'label': 'Z-score', 'shrink': 0.3, 'pad': 0.13})
    ax1.set_title('Motif Abundance per Sample (Z-score normalized)')
    # Force all ticks to appear
    ax1.set_xticks(np.arange(per_sample_abundance.T.shape[1]) + 0.5)
    ax1.set_xticklabels(per_sample_abundance.T.columns, rotation=90)

    # Add p-values to the right of the heatmap
    if motif_association_df is not None:
        for i in range(len(motif_association_df)):
            y_pos = i
            x_pos = len(per_sample_abundance) + 0.3
            p_val = motif_association_df.iloc[i]['p_value']
            ax1.text(x_pos, y_pos + 0.5, f"{p_val:.3f}", ha='left', va='center', fontsize=10,
                    transform=ax1.transData, color='red' if p_val < 0.05 else 'black')

    # Rename motifs to CN1, CN2, etc. based on actual motif IDs
    if prepare_motif_labels:
        motif_labels = [f"CN{motif+1}" for motif in per_sample_abundance.columns]
    else:
        motif_labels = per_sample_abundance.columns.tolist()
    ax1.set_yticklabels(motif_labels, rotation=0)

    # Rotate x-ticks and set alignment
    ax1.tick_params(axis='x', rotation=45)
    plt.setp(ax1.get_xticklabels(), ha='right')

    # Right figure: motif composition
    colors = ['white', '#FF0000']
    cmap = LinearSegmentedColormap.from_list('white_to_red', colors)
    cbar_label = 'Log count' if mode == 'log' else 'Raw'
    if mode == 'log':
        per_motif_composition = np.log2(1 + per_motif_composition)
    sns.heatmap(per_motif_composition, annot=False, vmin=0, cmap=cmap,
                xticklabels=1, yticklabels=1,
                cbar_kws={'label': cbar_label, 'shrink': 0.3, 'pad': 0.028}, ax=ax2)
    ax2.set_title(f'Motif Phenotype Composition',  pad=6)

    # Rename motifs to CN1, CN2, etc. on right heatmap
    ax2.set_yticklabels(motif_labels, rotation=0)

    # Rotate x-ticks on right heatmap and replace underscores with spaces
    ax2.tick_params(axis='x', rotation=45)
    xticklabels = ax2.get_xticklabels()
    ax2.set_xticklabels([label.get_text().replace('_', ' ') for label in xticklabels], rotation=45, ha='right')

    # Hide the empty dendrogram axes
    ax_dendro.axis('off')

    # Add overall figure title
    if treatment and m:
        fig.suptitle(f'Cellular Neighborhood Analysis - {treatment} Treatment (m={m})',
                    fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, top=0.96)
    plt.savefig(outfile, bbox_inches='tight', dpi=300)