
import argparse
from glob import glob
import os
import shutil
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations

from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

parser = argparse.ArgumentParser(description="Analyze permutation test results")
parser.add_argument('--selector', type=str, default='larc-samples/attraction_repulsion_results/*.csv')
parser.add_argument("--output-dir", type=str, default="attraction_repulsion_heatmaps_larc", help="Output directory for plots")
# Plotting options
parser.add_argument('--display-digits', action='store_true', help="Whether to display digits in heatmap cells")
args = parser.parse_args()

def draw_heatmap_v2(df):
    print("Drawing heatmaps...")
    print("Dataframe shape:", df.shape)
    print("Dataframe columns:", df.columns.tolist())
    phenotypes = sorted(df['source_phenotype'].unique())
    heatmap_res = {(a, b): [] for a, b in combinations(['Observed', 'Permuted'], 2)}
    for src_pheno in phenotypes:
        for dest_pheno in phenotypes:
            sub_df = df[(df['source_phenotype'] == src_pheno) & (df['target_phenotype'] == dest_pheno)]
            # Contrary to the script 45, not any groups is considered here
            observed_values = sub_df['observed_value'].values
            permuted_values = sub_df['permuted_average'].values
            for (g1, g1_df), (g2, g2_df) in [
                (("Observed", observed_values), ("Permuted", permuted_values))
            ]:
                # Do a t-test
                if len(g1_df) > 1 and len(g2_df) > 1:
                    t_stat, p_val = stats.ttest_ind(g1_df, g2_df, equal_var=False)
                    mean_1 = np.mean(g1_df)
                    mean_2 = np.mean(g2_df)
                    heatmap_res[(g1, g2)].append({
                        'source_phenotype': src_pheno,
                        'target_phenotype': dest_pheno,
                        't_statistic': t_stat,
                        'p_value': p_val,
                        'mean_1': mean_1,
                        'mean_2': mean_2
                    })

    heatmap_dfs = {}
    for (group1, group2), results_list in heatmap_res.items():
        if len(results_list) > 0:
            heatmap_dfs[(group1, group2)] = pd.DataFrame(results_list)

        # Create -log(p-value) heatmaps for each comparison
        for (group1, group2), results_df in heatmap_dfs.items():

            # Calculate conditional log(p-value) based on mean differences
            def calc_conditional_log_pvalue(row):
                p_val = row['p_value']
                mean_diff = row['mean_1'] - row['mean_2']
                log_p = np.log10(np.clip(p_val, 1e-300, None))
                if mean_diff > 0:
                    return log_p  # -log(p-value) for positive differences
                else:
                    return -log_p  # log(p-value) for negative differences

            results_df['conditional_log_pvalue'] = results_df.apply(calc_conditional_log_pvalue, axis=1)

            # Create pivot table for heatmap
            heatmap_data = results_df.pivot(index='source_phenotype',
                                            columns='target_phenotype',
                                            values='conditional_log_pvalue')

            # Fill NaN values with 0
            heatmap_data = heatmap_data.fillna(0)

            # Calculate max absolute value for symmetric colorscale
            max_abs_value = np.abs(heatmap_data.values).max()

            # Create custom diverging colormap
            custom_cmap = LinearSegmentedColormap.from_list('custom_diverging',
                                                             ['#c63461', '#ffffff', '#4287de'])

            # Pre-compute linkage matrix for synchronized clustering
            from scipy.cluster.hierarchy import linkage
            from scipy.spatial.distance import pdist

            # Since rows and columns have identical labels (phenotypes),
            # we can use the same linkage for both axes to ensure synchronized clustering
            distance_matrix = pdist(heatmap_data, metric='euclidean')
            linkage_matrix = linkage(distance_matrix, method='ward')

            # Create hierarchically clustered heatmap with synchronized clustering
            g = sns.clustermap(heatmap_data,
                               # Use the same linkage for both rows and columns
                               row_linkage=linkage_matrix,
                               col_linkage=linkage_matrix,

                               # Visual parameters
                               cmap=custom_cmap,  # Custom pink/white/blue colormap
                               center=0,  # Center white at zero
                               vmin=-max_abs_value,  # Symmetric minimum
                               vmax=max_abs_value,   # Symmetric maximum
                               annot=args.display_digits,  # Show/hide cell values based on flag
                               fmt='.2f' if args.display_digits else '',  # Two decimal places if showing

                               # Layout parameters
                               figsize=(12, 10),
                               row_cluster=True,  # Keep clustering
                               col_cluster=True,  # Keep clustering
                               linewidths=0.5,  # Grid lines between cells

                               # Color bar
                               cbar_kws={'label': 'Conditional log₁₀(p-value)'} if args.display_digits else {}
                               )

            # Hide dendrogram axes but keep clustering
            g.ax_row_dendrogram.set_visible(False)
            g.ax_col_dendrogram.set_visible(False)

            # Customize labels and title
            title = f'Heatmap of pairwise cell-cell interactions (blue) or avoidances (red)'
            g.fig.suptitle(title, y=0.95, fontsize=14, fontweight='bold')

            # Format tick labels
            plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=11.88)
            plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=11.88)

            # Replace underscores with spaces in tick labels
            x_labels = [label.get_text().replace('_', ' ') for label in g.ax_heatmap.get_xticklabels()]
            y_labels = [label.get_text().replace('_', ' ') for label in g.ax_heatmap.get_yticklabels()]
            g.ax_heatmap.set_xticklabels(x_labels)
            g.ax_heatmap.set_yticklabels(y_labels)

            # Add axis labels
            g.ax_heatmap.set_xlabel('Target Phenotype', fontsize=11, fontweight='bold')
            g.ax_heatmap.set_ylabel('Source Phenotype', fontsize=11, fontweight='bold')

            # Add "Interaction" and "Avoidance" labels to colorbar
            cbar = g.ax_cbar
            if not args.display_digits:
                cbar.set_yticks([])  # Hide colorbar ticks when not displaying digits
            cbar.text(0.5, 1.05, 'Interaction', ha='center', va='bottom',
                     transform=cbar.transAxes, fontsize=10, fontweight='bold')
            cbar.text(0.5, -0.05, 'Avoidance', ha='center', va='top',
                     transform=cbar.transAxes, fontsize=10, fontweight='bold')

            # Save the plot
            os.makedirs(args.output_dir, exist_ok=True)
            filename = f"{args.output_dir}/{group1}_vs_{group2}_pvalue_heatmap.png"
            g.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Saved clustered conditional log(p-value) heatmap: {filename}")


if __name__ == '__main__':
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    gathered_df = None

    # Check if selector contains spaces (multiple files) or is a glob pattern
    selector_stripped = args.selector.strip('"')
    if ' ' in selector_stripped:
        # Space-separated list of file paths
        file_list = selector_stripped.split(' ')
    else:
        # Glob pattern
        file_list = glob(selector_stripped)

    print(f"Processing {len(file_list)} files...")

    # Step 1. Gather all data in a dataframe
    for file_path in file_list:
        print(f"Reading file: {file_path}")
        # We don't consider the slug
        ar_df = pd.read_csv(file_path)
        for stat in ["mean"]: # Only the 'mean' statisitc is used
            subset = ar_df[ar_df['statistic'] == stat].copy()
            # using map_effect_sig is not used (contrary to the original script 45)
            gathered_df = pd.concat([gathered_df, subset]) if gathered_df is not None else subset

    # Step 2. Draw heatmaps for each statistic
    if gathered_df is not None:
        draw_heatmap_v2(gathered_df)
    else:
        print("ERROR: No data was gathered. Please check your input files.")
    print()
