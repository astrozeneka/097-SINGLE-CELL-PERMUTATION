import pandas as pd
import os
import argparse
from glob import glob

from utils import get_sample_slug_from_path, df2adata, class_to_pheno
import scimap as sm

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