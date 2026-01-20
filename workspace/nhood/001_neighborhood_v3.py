import argparse
import pandas as pd

from utils import df2adata
import scimap as sm
import os


parser = argparse.ArgumentParser(description="Filter neighborhood data based on classification sets")
parser.add_argument('--input', type=str, required=True, help="Merged CSV data file")
parser.add_argument('--output', type=str, required=True, help="Output CSV file")
parser.add_argument("--method", type=str, default="knn", choices=["knn", "radius"])
parser.add_argument('--knn-count', type=int, default=20)
parser.add_argument('--radius', type=float, default=50.0)
# Column selectors
parser.add_argument('--x-coord-column', type=str, default='Centroid X', help="Column name for X coordinate")
parser.add_argument('--y-coord-column', type=str, default='Centroid Y', help="Column name for Y coordinate")
parser.add_argument('--phenotype-column', type=str, default='Phenotype', help="Column name for Phenotype")
parser.add_argument('--imageid-column', type=str, default='imageid_', help="Column name for Image ID")
parser.add_argument('--region-column', type=str, default='Parent Region', help="Column name for Region")

args = parser.parse_args()

if __name__ == '__main__':
    df = pd.read_csv(args.input)

    # print provided options
    print(f"Provided options: " + " ".join([f"{k}={v}" for k, v in vars(args).items()]))

    # patch the coordinates
    df.rename({
        args.x_coord_column: 'Centroid X',
        args.y_coord_column: 'Centroid Y',
        # args.phenotype_column: args.phenotype_column,
        args.imageid_column: 'imageid'
    }, axis=1, inplace=True)

    all_adata = df2adata(
        df,
        imageid='imageid',
        additional_cols=['Phenotype'],
    )

    print("Defining neighborhood counts with k =", args.knn_count)
    if args.method == "radius":
        print("Using radius method with radius =", args.radius)
        all_adata = sm.tl.spatial_count(
            all_adata,
            phenotype=args.phenotype_column,
            method="radius",
            radius=args.radius
        )
    elif args.method == "knn":
        print("Using knn method with knn =", args.knn_count)
        all_adata = sm.tl.spatial_count(
            all_adata,
            phenotype=args.phenotype_column,
            method="knn",
            knn=args.knn_count
        )

    # Merge datas with same index (all_adata.uns['spatial_count'] and all_adata.obs[["Phenotype"]])
    spatial_count_df = all_adata.uns['spatial_count']
    phenotype_df = all_adata.obs[[args.phenotype_column, args.region_column, 'imageid', 'X_centroid', 'Y_centroid']]
    merged_df = spatial_count_df.join(phenotype_df)
    merged_df.rename(columns={'X_centroid': args.x_coord_column, 'Y_centroid': args.y_coord_column, 'imageid': args.imageid_column}, inplace=True)
    all_adata.uns['spatial_count'] = merged_df

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    all_adata.uns['spatial_count'].to_csv(args.output)
    print("Done, saved to", args.output)