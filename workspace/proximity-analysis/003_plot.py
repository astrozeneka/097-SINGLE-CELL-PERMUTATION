import argparse
import ast
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--tuples")
parser.add_argument("--group1_name", default="Group 1")
parser.add_argument("--group2_name", default="Group 2")
parser.add_argument("--output")
args = parser.parse_args()

pheno_tuples = [tuple(t.split("|")) for t in args.tuples.split(",")]

stats_df = pd.read_csv(args.input)

stacked = []
stats = []
for pheno_a, pheno_b in pheno_tuples:
    row = stats_df[(stats_df['Phenotype A'] == pheno_a) & (stats_df['Phenotype B'] == pheno_b)].iloc[0]
    group1_values = ast.literal_eval(row['Group1 Values'])
    group2_values = ast.literal_eval(row['Group2 Values'])
    for v in group1_values:
        stacked.append({'Group': args.group1_name, 'Phenotype A': pheno_a, 'Phenotype B': pheno_b, 'Interaction/mm²': v})
    for v in group2_values:
        stacked.append({'Group': args.group2_name, 'Phenotype A': pheno_a, 'Phenotype B': pheno_b, 'Interaction/mm²': v})
    stats.append({'Phenotype A': pheno_a, 'Phenotype B': pheno_b, 'p-value': row['p-value']})

stacked_df = pd.DataFrame(stacked)
stacked_df['Phenotype Tuple'] = stacked_df['Phenotype A'] + ' ↔ ' + stacked_df['Phenotype B']

fig, ax = plt.subplots(figsize=(2 * len(pheno_tuples), 3.8))
palette = {args.group1_name: '#80EF80', args.group2_name: '#E0115F'}
sns.boxplot(data=stacked_df, x='Phenotype Tuple', y='Interaction/mm²', hue='Group', ax=ax,
            palette=palette, showfliers=False, legend=False, width=0.8, gap=0.12)
sns.stripplot(data=stacked_df, x='Phenotype Tuple', y='Interaction/mm²', hue='Group',
              dodge=True, color='black', alpha=0.7, size=5, ax=ax, legend=False)

n_tuples = len(pheno_tuples)
y_pos = ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012
for i in range(n_tuples):
    ax.text(i - 0.2, y_pos, args.group1_name, ha='center', va='top', fontsize=9, transform=ax.transData)
    ax.text(i + 0.2, y_pos, args.group2_name, ha='center', va='top', fontsize=9, transform=ax.transData)

y_max = stacked_df['Interaction/mm²'].max()
y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
for idx, stat in enumerate(stats):
    x = idx
    y = y_max + y_range * 0.05
    ax.plot([x - 0.2, x - 0.2, x + 0.2, x + 0.2], [y, y + y_range * 0.02, y + y_range * 0.02, y], 'k-', linewidth=1)
    ax.text(x, y + y_range * 0.03, f"$p$={stat['p-value']:.3f}", ha='center', fontsize=8,
            color='red' if stat['p-value'] < 0.05 else 'black')

ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=8)
ax.tick_params(axis='x', pad=15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
ax.set_xlabel('')
plt.savefig(args.output)
