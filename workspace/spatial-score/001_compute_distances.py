import argparse
import os
import shutil
import zipfile
import pandas as pd
import anndata as ad
import scimap as sm

parser = argparse.ArgumentParser(description="Compute spatial distances, splitting by imageid")
parser.add_argument('--object-id', type=str, required=True)
parser.add_argument('--centroid-x', type=str, required=True)
parser.add_argument('--centroid-y', type=str, required=True)
parser.add_argument('--phenotype', type=str, required=True)
parser.add_argument('--imageid', type=str, required=True)
parser.add_argument('--output', type=str, default='data/spatial_distances_output.zip')
args, unknown = parser.parse_known_args()

# Collect all --inputN arguments
input_files = [v for k, v in zip(unknown[::2], unknown[1::2]) if k.startswith('--input')]


def df2adata(df, object_id_col, x_coord_col, y_coord_col, phenotype_col, imageid):
    expr = pd.DataFrame({'dummy': [0] * len(df)})
    meta = pd.DataFrame({
        'Object ID': df[object_id_col].values,
        'X_centroid': df[x_coord_col].values,
        'Y_centroid': df[y_coord_col].values,
        'Phenotype': df[phenotype_col].values,
        'imageid': imageid
    })
    adata = ad.AnnData(expr.values)
    adata.var.index = ['dummy']
    adata.obs = meta.set_index('Object ID')
    adata.raw = adata.copy()
    return adata


if __name__ == '__main__':
    tmp_dir = "spatial_distance_raw"
    os.makedirs(tmp_dir, exist_ok=True)

    for f in input_files:
        image_id = os.path.basename(f).replace('.csv', '')
        print(f"Processing {image_id} ({f})...")
        df = pd.read_csv(f)
        print(f"  {len(df)} cells")
        adata = df2adata(df, args.object_id, args.centroid_x, args.centroid_y, args.phenotype, image_id)

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
        output_df.to_csv(f"{tmp_dir}/{image_id}.csv", index=False)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with zipfile.ZipFile(args.output, 'w') as z:
        for f in os.listdir(tmp_dir):
            z.write(f"{tmp_dir}/{f}", f)
    shutil.rmtree(tmp_dir)
    print(f"Output saved to {args.output}")
    print("Done")
