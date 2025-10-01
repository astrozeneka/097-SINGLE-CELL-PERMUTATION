
import os
import argparse
from glob import glob
import shutil
import pandas as pd

parser = argparse.ArgumentParser(description="Analyze permutation test results")
parser.add_argument('--file-tuple', type=str, default="""[
    ["048a", "observed/spatial_distance_SAB_combined_zscore_048a.tif - 3 roi neighbourhoods.csv", "permuted/permutation_results_SAB_combined_zscore_048a.tif - 3 roi neighbourhoods.csv"],
    ["434r", "observed/spatial_distance_SAB_combined_zscore_434r 3 and 4 neighbourhood measurements.csv", "permuted/permutation_results_SAB_combined_zscore_434r 3 and 4 neighbourhood measurements.csv"],
    ["838r", "observed/spatial_distance_SAB_combined_zscore_838r.tif - 3 roi neighbourhood measurement.csv", "permuted/permutation_results_SAB_combined_zscore_838r.tif - 3 roi neighbourhood measurement.csv"],
    ["m043r", "observed/spatial_distance_SAB_combined_zscore_m043r.tif - 3 roi neighbourhood t1t7t9.csv", "permuted/permutation_results_SAB_combined_zscore_m043r.tif - 3 roi neighbourhood t1t7t9.csv"],
    ["m152p", "observed/spatial_distance_SAB_combined_zscore_m152p 3 and 8 roi neighbourhood measurement.csv", "permuted/permutation_results_SAB_combined_zscore_m152p 3 and 8 roi neighbourhood measurement.csv"],
    ["m573pl", "observed/spatial_distance_SAB_combined_zscore_m573pl.tif -3 roi neighbourhood.csv", "permuted/permutation_results_SAB_combined_zscore_m573pl.tif -3 roi neighbourhood.csv"]
]""")
parser.add_argument('--output-dir', type=str, default='attraction_repulsion_results')
args = parser.parse_args()

def get_file_slug_from_path(file_path):
    base_name = os.path.basename(file_path)
    if "permutation_results" in file_path:
        return base_name.replace("_spatial_distance_permutations.csv", "")
    if "spatial_distance_results" in file_path:
        return base_name.replace("_spatial_distances.csv", "")
    return None

# Prepare the file pair list
import json
file_pairs = json.loads(args.file_tuple)
file_pairs = {item[0]: (item[2], item[1]) for item in file_pairs}
#permutation_files = glob(args.permutation_selector.strip('"'))
#observed_files = glob(args.observed_selector.strip('"'))
#for pf in permutation_files:
#    slug = get_file_slug_from_path(pf)
    #    file_pairs[slug] = (pf, None)
#for of in observed_files:
#    slug = get_file_slug_from_path(of)
#    if slug in file_pairs:
#        file_pairs[slug] = (file_pairs[slug][0], of)

if __name__ == '__main__':

    # Delete dir if it exists
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    for sample, (permutation_file, observed_file) in file_pairs.items():
        print(f"Analyzing sample: {sample}")
        # Compute stats (means/median/p10/p90) for the observed file
        observed_data = pd.read_csv(observed_file, index_col=0)
    # Convert into a bi-indexed (src, dst) Dataframe, more efficient for computations
        # stat_df = pd.DataFrame(columns=["source_phenotype", "target_phenotype", "distance"])
        pheno_counts = observed_data['Phenotype'].value_counts()
        usable_phenotypes = pheno_counts[pheno_counts >= 10].index.tolist()
        usable_phenotypes = [u for u in usable_phenotypes if u != 'Other']

        observed_stats = []
        for src_pheno in usable_phenotypes:
            for dest_pheno in usable_phenotypes:
                sub_df = observed_data[[src_pheno, 'Phenotype']]
                sub_df = sub_df[sub_df['Phenotype'] == dest_pheno]
                stat_row = {
                    "source_phenotype": src_pheno,
                    "target_phenotype": dest_pheno,
                    "mean": sub_df[src_pheno].mean(),
                    "median": sub_df[src_pheno].median(),
                    "p10": sub_df[src_pheno].quantile(0.1),
                    "p90": sub_df[src_pheno].quantile(0.9),
                    "count": len(sub_df)
                }
                observed_stats.append(stat_row)
        observed_stats = pd.DataFrame(observed_stats)


        # Compute the p-value by counting the number of times the observed value is extreme than the randomized values
        sample_attraction_repulsion = []
        permutation_data_df = pd.read_csv(permutation_file)
        for src_pheno in usable_phenotypes:
            for dest_pheno in usable_phenotypes:
                sub_permu_df = permutation_data_df[
                    (permutation_data_df['source_phenotype'] == src_pheno) &
                    (permutation_data_df['target_phenotype'] == dest_pheno)
                ]
                for stat in ["mean", "median", "p10", "p90"]:
                    # For each statistics, compute empirical p-value
                    permu_avg = sub_permu_df[stat].mean()
                    curr_val = observed_stats[
                        (observed_stats['source_phenotype'] == src_pheno) &
                        (observed_stats['target_phenotype'] == dest_pheno)
                    ][stat].values[0] # This is expected to be a single value
                    # Count how many permuted situtation is more extreme than the observed value
                    if curr_val < permu_avg:
                        more_extreme_count = (sub_permu_df[stat] <= curr_val).sum()
                        p_val = more_extreme_count / len(sub_permu_df)
                        effect = "attraction"
                    else: # curr_val >= permu_avg
                        more_extreme_count = (sub_permu_df[stat] >= curr_val).sum()
                        p_val = more_extreme_count / len(sub_permu_df)
                        effect = "avoidance"
                    sample_attraction_repulsion.append({
                        "source_phenotype": src_pheno,
                        "target_phenotype": dest_pheno,
                        "statistic": stat,
                        "observed_value": curr_val,
                        "permuted_average": permu_avg,
                        "p_value": p_val,
                        "effect": effect
                    })

        sample_attraction_repulsion = pd.DataFrame(sample_attraction_repulsion)
        sample_attraction_repulsion.to_csv(f"attraction_repulsion_results/{sample}_attraction_repulsion.csv", index=False)