import pandas as pd
import os
import argparse
from glob import glob

from utils import class_to_pheno_v2, get_sample_slug_from_path, df2adata, class_to_pheno
import scimap as sm

#phenotypes = ['Tumor_cells', 'Tumor_proliferating_cells', 'B_cells', 'Vessels', 'Cytotoxic_T_cells', 'Helper_T_cells',
#              'Exhausted_Cytotoxic_T_cells', 'Exhausted_Helper_T_cells', 'Tregs', 'Other_T_cells',
#              'Antigen_presenting_cells', 'M1_Macrophages', 'M2_Macrophages', 'Other']

def process_single_file(file_path, n_permutations=25):
    df = pd.read_csv(file_path)
    df["Phenotype"] = df["Classification"].apply(class_to_pheno)
    # The total list of usable phenotype is all phenotype with at least 10 cells
    pheno_counts = df['Phenotype'].value_counts()
    usable_phenotypes = pheno_counts[pheno_counts >= 10].index.tolist()
    usable_phenotypes = [u for u in usable_phenotypes if u != 'Other']
    df = df[['Object ID', 'Parent Region', 'Parent Area µm^2', 'Parent Classification', 'Classification', 'Area', 'Centroid X',
             'Centroid Y', 'Phenotype']]

    sample_slug = get_sample_slug_from_path(file_path)
    # Subsample for quick testing
    # df = df.head(5000)  # TODO, remove this line for full analysis
    all_stat = None
    for i in range(n_permutations):
        print("Processing permutation", i + 1, "of", n_permutations, end='\r')
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

    os.makedirs(f"permutation_results/n{n_permutations}", exist_ok=True)
    all_stat.to_csv(f"permutation_results/n{n_permutations}/{sample_slug}_spatial_distance_permutations.csv", index=False)
    print()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute spatial interaction metrics between phenotypes")
    parser.add_argument('--selector', type=str,
                        default='../013-LARC/*/proj-qupath-alignment-*_5Cycles_aligned_tschnm+class/cell_data_with_regions_and_parent_areas*.csv')
    parser.add_argument("--slug", type=str, default="_T", help="Sample slug to process")
    parser.add_argument("--n_permutations", type=int, default=25, help="Number of permutations to perform")
    args = parser.parse_args()

    file_list = glob(args.selector.strip('"'))
    file_list = [f for f in file_list if args.slug in f]

    for file_path in file_list:
        print(f"Processing file: {file_path}")
        process_single_file(file_path)