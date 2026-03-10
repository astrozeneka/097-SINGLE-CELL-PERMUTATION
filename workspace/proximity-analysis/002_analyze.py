import argparse
import pandas as pd
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument("--group1_inputs")
parser.add_argument("--group2_inputs")
parser.add_argument("--output")
args = parser.parse_args()

group1_paths = args.group1_inputs.split(",") if args.group1_inputs else []
group2_paths = args.group2_inputs.split(",") if args.group2_inputs else []

def load_matrix(path):
    return pd.read_csv(path, index_col=0)

group1_dfs = [load_matrix(p) for p in group1_paths]
group2_dfs = [load_matrix(p) for p in group2_paths]

pheno_list = sorted(set().union(*[set(df.index) for df in group1_dfs + group2_dfs]))

def reindex(df):
    return df.reindex(index=pheno_list, columns=pheno_list, fill_value=0)

group1_dfs = [reindex(df) for df in group1_dfs]
group2_dfs = [reindex(df) for df in group2_dfs]

stats = []
for pheno_a in pheno_list:
    for pheno_b in pheno_list:
        if pheno_list.index(pheno_b) < pheno_list.index(pheno_a):
            continue
        group1_values = [float(df.loc[pheno_a, pheno_b]) for df in group1_dfs]
        group2_values = [float(df.loc[pheno_a, pheno_b]) for df in group2_dfs]
        mean_group1 = sum(group1_values) / len(group1_values) if group1_values else 0
        mean_group2 = sum(group2_values) / len(group2_values) if group2_values else 0
        u, p = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
        stats.append({
            'Phenotype A': pheno_a,
            'Phenotype B': pheno_b,
            'Mean Group 1': mean_group1,
            'Mean Group 2': mean_group2,
            'Group1 Values': group1_values,
            'Group2 Values': group2_values,
            'Mann-Whitney U': u,
            'p-value': p,
        })

# sort by ascending p
stats.sort(key=lambda x: x['p-value'])
pd.DataFrame(stats).to_csv(args.output, index=False)
