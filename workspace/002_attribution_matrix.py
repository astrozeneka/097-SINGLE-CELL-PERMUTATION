

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore

parser = argparse.ArgumentParser("Build the composition matrix for one sample")
parser.add_argument('--input', type=str, required=True, help="Input CSV file with classified cells")
parser.add_argument('--output', type=str, required=True, help="Output PNG file for the heatmap")

# Optional arguments can be added here
parser.add_argument('--phenotypes', type=str, default='Phenotype')
parser.add_argument('--phenotype', type=str, default='Phenotype')
parser.add_argument('--markers', type=str)

# CD68,CD163,CD206,CD11c,PDL1,CD3,CD4,FOXP3,CD8,PD1,CD20,CD31,Ki67
# Mac,Mac_CD163,Mac_CD206,Mac_CD163_CD206,DC,T,T_CD4,T_CD8,Treg,T_PD1,T_CD4_PD1,T_CD8_PD1,B,Endothelial cells
# TAM,DC,T,Helper T,T PD1,Tregs,Endothelial cells,Ki67 proliferating cells

args = parser.parse_args()

def plot_heatmap(phenotype_means, save_path, phenotype_list):
    """Plot a single heatmap of phenotype means"""
    # Reorder phenotypes according to phenotype_list
    available_phenotypes = [p for p in phenotype_list if p in phenotype_means.index]
    phenotype_means_reordered = phenotype_means.reindex(index=available_phenotypes)

    # Calculate z-scores if we have enough phenotypes
    if len(phenotype_means_reordered) > 1:
        phenotype_zscore = phenotype_means_reordered.apply(lambda x: zscore(x, nan_policy='omit'), axis=0)
        cbar_label = 'Mean Intensity (Z-score)'
    else:
        phenotype_zscore = phenotype_means_reordered
        cbar_label = 'Mean Intensity'

    # Create heatmap
    plt.figure(figsize=(10, 8))

    # Calculate symmetric vmin/vmax for colorbar
    max_abs_value = max(abs(phenotype_zscore.min().min()), abs(phenotype_zscore.max().max()))

    ax = sns.heatmap(phenotype_zscore, cmap='bwr', center=0, annot=False,
                vmin=-max_abs_value, vmax=max_abs_value,
                cbar_kws={'label': cbar_label, 'shrink': 0.18, 'aspect': 4})

    # Add border to colorbar
    cbar = ax.collections[0].colorbar
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(1)

    plt.title('')
    plt.xlabel('Markers')
    plt.ylabel('Phenotypes')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    df = pd.read_csv(args.input, index_col=0)
    marker_list = args.markers.split(',')
    marker_list = [a.strip() for a in marker_list]

    phenotype_list = args.phenotypes.split(',')
    phenotype_list = [a.strip() for a in phenotype_list]

    # Calculate aggregated mean: average intensity by phenotype
    phenotype_means = df.groupby(args.phenotype)[marker_list].mean()

    print("Aggregated mean staining intensity by phenotype:")
    print(phenotype_means)
    print()

    # Save to CSV for inspection
    output_csv = args.output.replace('.png', '_means.csv')
    phenotype_means.to_csv(output_csv)
    print(f"Saved aggregated means to: {output_csv}")

    # Plot heatmap
    plot_heatmap(phenotype_means, args.output, phenotype_list)
    print(f"Saved heatmap to: {args.output}")