
import argparse

import pandas as pd
from sklearn.cluster import KMeans

parser = argparse.ArgumentParser(description="Cluster neighborhood into motifs.")
parser.add_argument('--input', type=str, default="nhood_results/all_k12.csv", help="Neighborhood results CSV file")
parser.add_argument('--n-motifs', type=int, default=24, help="Number of motifs to cluster")
parser.add_argument('--knn-seed', type=int, default=42, help="Random seed for KMeans")
parser.add_argument("--output", type=str, required=True, help="Output file name")
args = parser.parse_args()

if __name__ == '__main__':
    nhood_res = pd.read_csv(args.input, index_col=0)
    nhood_vec = nhood_res.drop(columns=['Phenotype', 'imageid'])
    kmeans = KMeans(n_clusters=args.n_motifs, random_state=args.knn_seed)
    motifs = kmeans.fit_predict(nhood_vec)
    nhood_res['motif'] = motifs
    # save to output
    nhood_res.to_csv(args.output)