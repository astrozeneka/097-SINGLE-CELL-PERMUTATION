import argparse
import os
import shutil
import zipfile
from scipy.spatial import KDTree
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("--inputs")
parser.add_argument("--centroid_x_col")
parser.add_argument("--centroid_y_col")
parser.add_argument("--parent_area_col")
parser.add_argument("--parent_region_col")
parser.add_argument("--cell_type_col")
parser.add_argument("--pixel_size", type=float, default=0.5)
parser.add_argument("--output")
args = parser.parse_args()

inputs = args.inputs.split(",") if args.inputs else []


def compute_cell_cell_interactions(df, sample_name, radius_in_micrometer=15, pixel_size=0.5, x_col='Centroid X', y_col='Centroid Y', phenotype_col='Phenotype', parent_region_col='Parent Region', parent_area_col='Parent Area µm^2'):
    radius_in_pixels = radius_in_micrometer / pixel_size
    pheno_list = df[phenotype_col].unique()
    print(f"Processing sample: {sample_name}")
    parent_regions = df[[parent_region_col, parent_area_col]].drop_duplicates()
    total_area_um2 = parent_regions[parent_area_col].sum()
    total_area_mm2 = total_area_um2 / 1e6

    coords = df[[x_col, y_col]].values
    phenotypes = df[phenotype_col].values
    tree = KDTree(coords)
    pairs = tree.query_pairs(r=radius_in_pixels, output_type='ndarray')
    matrix = pd.DataFrame(0, index=pheno_list, columns=pheno_list)
    if len(pairs) > 0:
        phenos_i = phenotypes[pairs[:, 0]]
        phenos_j = phenotypes[pairs[:, 1]]
        pair_count = Counter(zip(phenos_i, phenos_j))
        for (pheno_a, pheno_b), count in pair_count.items():
            matrix.loc[pheno_a, pheno_b] += count
            if pheno_a != pheno_b:
                matrix.loc[pheno_b, pheno_a] += count
    matrix = matrix / total_area_mm2
    os.makedirs("interaction_matrix_raw", exist_ok=True)
    matrix.to_csv(f"interaction_matrix_raw/{sample_name}.csv")


for f in inputs:
    sample_name = '_'.join(os.path.splitext(os.path.basename(f))[0].split('_')[:-1])
    df = pd.read_csv(f)
    compute_cell_cell_interactions(
        df,
        sample_name=sample_name,
        x_col=args.centroid_x_col,
        y_col=args.centroid_y_col,
        parent_area_col=args.parent_area_col,
        parent_region_col=args.parent_region_col,
        phenotype_col=args.cell_type_col,
        pixel_size=args.pixel_size,
    )

if args.output:
    output_path = args.output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with zipfile.ZipFile(output_path, 'w') as z:
        for f in os.listdir("interaction_matrix_raw"):
            z.write(f"interaction_matrix_raw/{f}", f)
    shutil.rmtree("interaction_matrix_raw")
    print(f"Output saved to {output_path}")
