import argparse
import scimap as sm
import os
from glob import glob
import pandas as pd
import anndata as ad

def df2adata(df, imageid="sample", additional_cols=[]):

    feature_cols = [
        "DAPI1",
        "Ki67",
        "-",
        "CD11c",
        "DAPI2",
        "CD3d",
        "--",
        "MHCII",
        "DAPI3",
        "CD68",
        "CD8",
        "PD1",
        "DAPI4",
        "FOXP3",
        "CD4",
        "CD20",
        "DAPI5",
        "PanCK",
        "CD163",
        "CD31"
    ]
    for col in feature_cols:
        if col + ': Cell: Mean' in df.columns:
            # rename
            df.rename(columns={col + ': Cell: Mean': col}, inplace=True)

    # Define metadata columns - mapping to your new data structure
    meta_cols = [
        'Object ID',  # CellID
        'Centroid X',  # X_centroid
        'Centroid Y',  # Y_centroid
        'Parent Region',  # Additional metadata
        'Parent Area µm^2', # Parent Area
        'Parent Classification',
        'Classification',
        'Area',
        'imageid'  # Sample slug
    ] + additional_cols

    # Check which feature columns are actually present in the data
    available_features = [col for col in feature_cols if col in df.columns]
    missing_features = [col for col in feature_cols if col not in df.columns]

    print(f"Available features: {available_features}")
    #if missing_features:
    #    print(f"Missing features: {missing_features}")

    # Check metadata columns
    available_meta = [col for col in meta_cols if col in df.columns]
    missing_meta = [col for col in meta_cols if col not in df.columns]

    print(f"Available metadata columns: {available_meta}")
    #if missing_meta:
    #    print(f"Missing metadata: {missing_meta}")

    # Extract expression data and metadata
    expr = df[available_features].copy()
    meta = df[available_meta].copy()

    # Create simplified marker names for AnnData (remove ": Cell: Mean" suffix)
    marker_names = [col.replace(': Cell: Mean', '') for col in available_features]

    # Create AnnData object
    adata = ad.AnnData(expr.values)
    adata.var.index = marker_names  # simplified marker names

    # Set up observations (cells) metadata
    # Use Object ID as the index
    meta_indexed = meta.set_index('Object ID')
    adata.obs = meta_indexed

    # Rename coordinate columns to match scimap expectations
    if 'Centroid X' in adata.obs.columns:
        adata.obs['X_centroid'] = adata.obs['Centroid X']
    if 'Centroid Y' in adata.obs.columns:
        adata.obs['Y_centroid'] = adata.obs['Centroid Y']

    # Add sample slug to observations
    adata.obs['imageid'] = adata.obs[imageid]

    # Add all markers to uns
    adata.uns['all_markers'] = list(adata.var.index)

    # Add copy of the original data to raw
    adata.raw = adata.copy()

    return adata

parser = argparse.ArgumentParser(description="Calculate spatial distances between phenotypes")
parser.add_argument('--input', type=str, required=True, help='Input CSV file path')
parser.add_argument('--output', type=str, required=True, help='Output CSV file path')
parser.add_argument('--slug', type=str, required=True, help='Sample slug to assign as imageid')
# Additional optional parameters
parser.add_argument('--phenotype', type=str, default='Phenotype')
parser.add_argument('--x-coord', type=str, default='Centroid X')
parser.add_argument('--y-coord', type=str, default='Centroid Y')
args = parser.parse_args()

if __name__ == '__main__':
    print("Computing spatial distances for real tissue data")
    df = pd.read_csv(args.input)
    # Rename phenotype column if necessary
    if args.phenotype not in df.columns:
        raise ValueError(f"Phenotype column '{args.phenotype}' not found in data")
    df = df.rename(columns={args.phenotype: 'Phenotype'})
    df['imageid'] = args.slug

    adata = df2adata(df, imageid='imageid', additional_cols=['Phenotype'])
    adata = sm.tl.spatial_distance(
        adata,
        x_coordinate=args.x_coord,
        y_coordinate=args.y_coord,
        phenotype='Phenotype',
        imageid='imageid'
    )

    # Prepare and save csv
    adata.uns['spatial_distance'].to_csv(args.output)
    print("Output data saved to:", args.output)
    print("Done")
