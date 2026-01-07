
import argparse
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from plot_utils import create_combined_plot_v2_without_group

parser = argparse.ArgumentParser(description="Visualize neighborhood content")
parser.add_argument('--input', type=str, required=True, help="Input stacked neighborhood interaction CSV file including the motif information")
parser.add_argument("--output", type=str, required=True, help="Output figure file name")
parser.add_argument("--fig-width", type=float, default=12, help="Figure width")
args = parser.parse_args()

if __name__ == '__main__':
    nhood_res = pd.read_csv(args.input, index_col=0)
    print()

    # Code reusage/even copy-paste can help to increase coding efficiency
    clustered = nhood_res.copy()

    # A1. Motif abundance (ratio) per sample (aggregate + average + z-score)
    per_sample_abundance = clustered.groupby(['imageid', 'motif']).size().unstack(fill_value=0)
    per_sample_abundance = per_sample_abundance.div(per_sample_abundance.sum(axis=1), axis=0)
    per_sample_abundance = (per_sample_abundance - per_sample_abundance.mean()) / per_sample_abundance.std()

    # No need to organize
    # Might be included in the v2

    # B. Not implemented
    # Might be included in the v2

    # C1. Motif phenotype comosition
    per_motif_composition = clustered.groupby(['motif', 'Phenotype']).size().unstack(fill_value=0)
    per_motif_composition = per_motif_composition.div(per_motif_composition.sum(axis=1), axis=0)
    # Cluster rows and columns
    row_linkage = linkage(pdist(per_motif_composition, metric='euclidean'), method='ward')
    row_order = dendrogram(row_linkage, no_plot=True)['leaves']
    col_linkage = linkage(pdist(per_motif_composition.T, metric='euclidean'), method='ward')
    col_order = dendrogram(col_linkage, no_plot=True)['leaves']
    per_motif_composition = per_motif_composition.iloc[row_order, col_order]

    # reorder the per_sample_abundance columns to match the motif order in per_motif_composition
    per_sample_abundance = per_sample_abundance[per_motif_composition.index]

    # Call the function without response/group information
    create_combined_plot_v2_without_group(
        per_sample_abundance=per_sample_abundance,
        motif_association_df=None,
        per_motif_composition=per_motif_composition,
        outfile=args.output,
        treatment='All',
        m=None,
        mode='raw',
        prepare_motif_labels=True,
        fig_width=args.fig_width
    )