import argparse
import pandas as pd
import anndata as ad
import scimap as sm

parser = argparse.ArgumentParser(description="Compute spatial distances, splitting by imageid")
parser.add_argument('--object-id', type=str, required=True)
parser.add_argument('--centroid-x', type=str, required=True)
parser.add_argument('--centroid-y', type=str, required=True)
parser.add_argument('--phenotype', type=str, required=True)
parser.add_argument('--imageid', type=str, required=True)
parser.add_argument('--output', type=str, default='data/spatial_distances_output.csv')
args, unknown = parser.parse_known_args()

# Collect all --inputN arguments
input_files = [v for k, v in zip(unknown[::2], unknown[1::2]) if k.startswith('--input')]


def df2adata(df, object_id_col, x_coord_col, y_coord_col, phenotype_col, imageid_col='imageid'):
    expr = pd.DataFrame({'dummy': [0] * len(df)})
    meta = pd.DataFrame({
        'Object ID': df[object_id_col].values,
        'X_centroid': df[x_coord_col].values,
        'Y_centroid': df[y_coord_col].values,
        'Phenotype': df[phenotype_col].values,
        'imageid': df[imageid_col].values
    })
    adata = ad.AnnData(expr.values)
    adata.var.index = ['dummy']
    adata.obs = meta.set_index('Object ID')
    adata.raw = adata.copy()
    return adata


if __name__ == '__main__':
    # Merge all input files
    dfs = []
    for f in input_files:
        print(f"Reading {f}...")
        dfs.append(pd.read_csv(f))
    df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(df)} cells total")

    results = []
    for image_id, group_df in df.groupby(args.imageid):
        print(f"Processing imageid={image_id} ({len(group_df)} cells)...")
        adata = df2adata(group_df, args.object_id, args.centroid_x, args.centroid_y, args.phenotype, args.imageid)

        adata = sm.tl.spatial_distance(
            adata,
            x_coordinate='X_centroid',
            y_coordinate='Y_centroid',
            phenotype='Phenotype',
            imageid='imageid'
        )

        output_df = adata.uns['spatial_distance'].copy()
        output_df = output_df.assign(Phenotype=adata.obs['Phenotype'].values)
        output_df['imageid'] = image_id
        results.append(output_df)

    merged = pd.concat(results, ignore_index=True)
    merged.to_csv(args.output, index=False)
    print(f"Saved {len(merged)} rows to {args.output}")
    print("Done")
