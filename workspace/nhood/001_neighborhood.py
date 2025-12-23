import argparse
import pandas as pd

from utils import df2adata
# import scimap as sm
import os


parser = argparse.ArgumentParser(description="Filter neighborhood data based on classification sets")
parser.add_argument('--input', type=str, required=True, help="Merged CSV data file")
parser.add_argument('--output', type=str, required=True, help="Output CSV file")
parser.add_argument('--knn-count', type=int, default=20)
args = parser.parse_args()

if __name__ == '__main__':
    df = pd.read_csv(args.input)
    
    all_adata = df2adata(
        df,
        imageid="imageid",
        additional_cols=["Phenotype"]
    )

    print("Defining neighborhood counts with k =", args.knn_count)
    all_adata = sm.tl.spatial_count(
        all_adata,
        phenotype="Phenotype",
        method="knn",
        knn=args.knn_count
    )

    # Merge datas with same index (all_adata.uns['spatial_count'] and all_adata.obs[["Phenotype"]])
    spatial_count_df = all_adata.uns['spatial_count']
    phenotype_df = all_adata.obs[["Phenotype", "imageid"]]
    merged_df = spatial_count_df.join(phenotype_df)
    all_adata.uns['spatial_count'] = merged_df

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    all_adata.uns['spatial_count'].to_csv(args.output)
    print("Done, saved to", args.output)