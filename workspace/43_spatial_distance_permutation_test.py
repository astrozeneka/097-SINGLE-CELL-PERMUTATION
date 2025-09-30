import pandas as pd
import os
import argparse
from glob import glob

from utils import get_sample_slug_from_path, class_to_pheno
import scimap as sm
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

#phenotypes = ['Tumor_cells', 'Tumor_proliferating_cells', 'B_cells', 'Vessels', 'Cytotoxic_T_cells', 'Helper_T_cells',
#              'Exhausted_Cytotoxic_T_cells', 'Exhausted_Helper_T_cells', 'Tregs', 'Other_T_cells',
#              'Antigen_presenting_cells', 'M1_Macrophages', 'M2_Macrophages', 'Other']

def process_single_file(file_path, output_path, sample_slug, n_permutations=25):
    df = pd.read_csv(file_path)
    # The total list of usable phenotype is all phenotype with at least 10 cells
    pheno_counts = df['Phenotype'].value_counts()
    usable_phenotypes = pheno_counts[pheno_counts >= 10].index.tolist()
    usable_phenotypes = [u for u in usable_phenotypes if u != 'Other']
    df = df[['Object ID', 'Parent Region', 'Parent Area µm^2', 'Centroid X',
             'Centroid Y', 'Phenotype']]

    # Subsample for quick testing
    # df = df.head(5000)  # TODO, remove this line for full analysis
    all_stat = None
    for i in range(n_permutations):
        print("")
        print("** Processing permutation", i + 1, "of", n_permutations)
        df_permuted = df.copy()
        df_permuted['Phenotype'] = df_permuted['Phenotype'].sample(frac=1).values
        df_permuted['imageid'] = sample_slug
        adata = df2adata(df_permuted, imageid='imageid', additional_cols=['Phenotype'])
        # Calculate spatial distances
        adata = sm.tl.spatial_distance(
            adata,
            x_coordinate='X_centroid',
            y_coordinate='Y_centroid',
            phenotype='Phenotype',
            imageid='imageid'
        )

        # Step 2. Get statistics from each cell src/target pair (mean, median, p10, p90)
        spatial_distance_df = adata.uns['spatial_distance']
        spatial_distance_df = spatial_distance_df.assign(
            Phenotype=adata.obs['Phenotype']
        )
        for src_pheno in usable_phenotypes:
            for dest_pheno in usable_phenotypes:
                sub_df = spatial_distance_df[[src_pheno, 'Phenotype']]
                sub_df = sub_df[sub_df['Phenotype'] == dest_pheno]
                stat_row = {
                    'permutation': i,
                    'source_phenotype': src_pheno,
                    'target_phenotype': dest_pheno,
                    'mean': sub_df[src_pheno].mean(),
                    'median': sub_df[src_pheno].median(),
                    'p10': sub_df[src_pheno].quantile(0.1),
                    'p90': sub_df[src_pheno].quantile(0.9),
                    'count': len(sub_df),
                    'Sample': sample_slug
                }
                if all_stat is None:
                    all_stat = pd.DataFrame([stat_row])
                else:
                    all_stat = pd.concat([all_stat, pd.DataFrame([stat_row])], ignore_index=True)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    all_stat.to_csv(output_path, index=False)
    print()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute spatial interaction metrics between phenotypes")
    parser.add_argument('--input', type=str, required=True, help='Input CSV file path')
    parser.add_argument("--sample", type=str, default=None, help="Sample ID")
    parser.add_argument('--output', type=str, required=True, help='Output CSV file path')
    parser.add_argument("--n_permutations", type=int, default=25, help="Number of permutations to perform")
    args = parser.parse_args()

    print(f"Processing file: {args.input}")
    process_single_file(args.input, args.output, args.sample, args.n_permutations)