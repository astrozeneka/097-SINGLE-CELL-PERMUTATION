import pandas as pd
import anndata as ad
import os

phenotype_df = pd.read_csv("phenotypes.csv")
phenotype_list = phenotype_df['Phenotype'].tolist()
phenotypes = list(set(phenotype_list))

def class_to_pheno(classification):
    phenotype = phenotype_df[['Phenotype', 'Classification']].loc[phenotype_df['Classification'] == classification]
    if not phenotype.empty:
        return phenotype['Phenotype'].values[0]
    else:
        return 'Other'

_phenotype_dict = None
def class_to_pheno_v2(classification, phenotype_csv="../080-QUPATH_VS_SCIMAP/phenotypes.csv"):
    global _phenotype_dict
    if _phenotype_dict is None:
        df = pd.read_csv(phenotype_csv)
        _phenotype_dict = pd.Series(df['Phenotype'].values, index=df['Classification']).to_dict()
    return _phenotype_dict.get(classification, 'Other')


def load_patients_df(data_path="patients.csv"):
    df = pd.read_csv(data_path)
    return df

def get_sample_slug_from_path(file_path):
    if "proj-qupath-alignment-" in file_path:  # If the data source is in the desktop pipeline
        sample_slug = file_path[file_path.find('proj-qupath-alignment-'):].split("+")[0] \
            .replace("proj-qupath-alignment-", "").replace("_5Cycles_aligned_tschnm", "")
    else:  # in notebook pipeline
        sample_slug = os.path.basename(file_path).replace('cell_data_with_regions_and_parent_areas_', '') \
            .replace('.csv', '')
        if "_base" in sample_slug or "_T" in sample_slug:
            sample_slug = sample_slug.split("_base")[0].split("_T")[0]
    return sample_slug

def df2adata(df, imageid="sample", additional_cols=[]):

    print(f"Loaded data with shape: {df.shape}")
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

    print(f"Available feature columns: {len(available_features)}")
    print(f"Available features: {available_features}")
    if missing_features:
        print(f"Missing features: {missing_features}")

    # Check metadata columns
    available_meta = [col for col in meta_cols if col in df.columns]
    missing_meta = [col for col in meta_cols if col not in df.columns]

    print(f"Available metadata columns: {available_meta}")
    if missing_meta:
        print(f"Missing metadata: {missing_meta}")

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

def get_round_from_slug(slug):
    items = {'S19_12126B1': 'Round1', 'S19_25142A1': 'Round1', 'S19_31776B1': 'Round1', 'S20_12903C1': 'Round1', 'S20_19504A1': 'Round3', 'S20_2317A9': 'Round3', 'S20_23887A1': 'Round1', 'S20_25858B1': 'Round3', 'S20_8961A8': 'Round1', 'S21_12126A1': 'Round2', 'S21_12210A3': 'Round2', 'S21_12602A10': 'Round2', 'S21_1580B1': 'Round2', 'S21_18334A9': 'Round1', 'S21_19043A1': 'Round3', 'S21_8615C1': 'Round2', 'S21_8643A13': 'Round1', 'S22_10684A1': 'Round2', 'S22_12362A1': 'Round1', 'S22_20017A4': 'Round3', 'S22_7776A7': 'Round2', 'S23_2629A6': 'Round1', 'S23_3055A18': 'Round2', 'SP18_2074B1': 'Round3', 'SP19_4321C1': 'Round1', 'SP21_204B1': 'Round1', 'SP22_1201A1': 'Round2', 'SP22_370B1': 'Round2'}
    return items[slug]