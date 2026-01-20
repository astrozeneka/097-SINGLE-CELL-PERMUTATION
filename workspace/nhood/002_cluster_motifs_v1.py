
import argparse

import pandas as pd
from sklearn.cluster import KMeans

parser = argparse.ArgumentParser(description="Cluster neighborhood into motifs.")
parser.add_argument('--input', type=str, default="nhood_results/all_k12.csv", help="Neighborhood results CSV file")
parser.add_argument('--n-motifs', type=int, default=24, help="Number of motifs to cluster")
parser.add_argument('--knn-seed', type=int, default=42, help="Random seed for KMeans")
parser.add_argument("--output", type=str, required=True, help="Output file name")

# Column selectors
parser.add_argument('--phenotype-column', type=str, default='Phenotype', help="Column name for Phenotype")
parser.add_argument('--imageid-column', type=str, default='imageid', help="Column name for Image ID")
parser.add_argument('--x-coord-column', type=str, default='Centroid X', help="Column name for X coordinate")
parser.add_argument('--y-coord-column', type=str, default='Centroid Y', help="Column name for Y coordinate")
parser.add_argument('--region-column', type=str, default='Parent Region', help="Column name for Region")

args = parser.parse_args()

if __name__ == '__main__':
    nhood_res = pd.read_csv(args.input, index_col=0)

    # Columns to exclude from clustering (informative/metadata columns)
    exclude_cols = [args.phenotype_column, args.imageid_column,
                    args.x_coord_column, args.y_coord_column,
                    args.region_column]

    nhood_vec = nhood_res.drop(columns=exclude_cols, errors='ignore')
    kmeans = KMeans(n_clusters=args.n_motifs, random_state=args.knn_seed)
    motifs = kmeans.fit_predict(nhood_vec)
    nhood_res['motif'] = motifs
    # save to output
    nhood_res.to_csv(args.output)