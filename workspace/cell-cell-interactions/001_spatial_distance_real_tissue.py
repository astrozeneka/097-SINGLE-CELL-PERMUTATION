import argparse
import scimap as sm
import pandas as pd
import anndata as ad


def df2adata(df, object_id_col, x_coord_col, y_coord_col, phenotype_col, imageid_col='imageid'):
    """
    Convert a DataFrame to AnnData format suitable for scimap spatial analysis.

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with cell data
    object_id_col : str
        Column name for cell/object IDs
    x_coord_col : str
        Column name for X coordinates
    y_coord_col : str
        Column name for Y coordinates
    phenotype_col : str
        Column name for phenotype/cell type
    imageid_col : str
        Column name for image/sample ID (default: 'imageid')

    Returns:
    --------
    ad.AnnData
        AnnData object ready for scimap spatial_distance
    """
    # Create minimal expression matrix (required by AnnData but not used for spatial distance)
    # Use a single dummy feature column
    expr = pd.DataFrame({'dummy': [0] * len(df)})

    # Create metadata DataFrame
    meta = pd.DataFrame({
        'Object ID': df[object_id_col].values,
        'X_centroid': df[x_coord_col].values,
        'Y_centroid': df[y_coord_col].values,
        'Phenotype': df[phenotype_col].values,
        'imageid': df[imageid_col].values
    })

    # Create AnnData object
    adata = ad.AnnData(expr.values)
    adata.var.index = ['dummy']

    # Set up observations (cells) metadata with Object ID as index
    meta_indexed = meta.set_index('Object ID')
    adata.obs = meta_indexed

    # Add copy of the original data to raw
    adata.raw = adata.copy()

    return adata


parser = argparse.ArgumentParser(description="Calculate spatial distances between phenotypes")
parser.add_argument('--input', type=str, required=True, help='Input CSV file path')
parser.add_argument('--output', type=str, required=True, help='Output CSV file path')
parser.add_argument('--object-id-column', type=str, required=True, help='Column name for cell/object IDs')
parser.add_argument('--x-coord-column', type=str, required=True, help='Column name for X coordinates')
parser.add_argument('--y-coord-column', type=str, required=True, help='Column name for Y coordinates')
parser.add_argument('--phenotype-column', type=str, required=True, help='Column name for phenotype/cell type')
parser.add_argument('--sample-name', type=str, required=True, help='Sample name to assign as imageid')
args = parser.parse_args()

if __name__ == '__main__':
    print("Computing spatial distances for real tissue data")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")

    # Read input CSV
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} cells")

    # Validate required columns exist
    required_cols = [args.object_id_column, args.x_coord_column, args.y_coord_column, args.phenotype_column]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in input data: {missing_cols}")

    # Add imageid column with sample name
    df['imageid'] = args.sample_name

    # Convert to AnnData
    adata = df2adata(
        df,
        object_id_col=args.object_id_column,
        x_coord_col=args.x_coord_column,
        y_coord_col=args.y_coord_column,
        phenotype_col=args.phenotype_column,
        imageid_col='imageid'
    )

    print(f"Computing spatial distances for {len(adata.obs['Phenotype'].unique())} phenotypes...")

    # Calculate spatial distances using scimap
    adata = sm.tl.spatial_distance(
        adata,
        x_coordinate='X_centroid',
        y_coordinate='Y_centroid',
        phenotype='Phenotype',
        imageid='imageid'
    )

    # Prepare output: add phenotype column back to spatial_distance result
    output_df = adata.uns['spatial_distance'].copy()
    output_df = output_df.assign(Phenotype=adata.obs['Phenotype'])

    # Rename phenotype column back to match input column name
    output_df = output_df.rename(columns={'Phenotype': args.phenotype_column})

    # Save to CSV
    output_df.to_csv(args.output)
    print(f"Output data saved to: {args.output}")
    print("Done")
