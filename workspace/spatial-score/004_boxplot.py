import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--tuples', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

df = pd.read_csv(args.input)
group_names = [c[:-5] for c in df.columns if c.endswith('_vals')]

selected = [t.split('|') for t in args.tuples.split(',')]
df = df[df.apply(lambda r: [r['target'], r['effector']] in selected, axis=1)].reset_index(drop=True)

rows = []
for _, row in df.iterrows():
    label = f"{row['target']} → {row['effector']}"
    for g in group_names:
        vals = [float(v) for v in str(row[f'{g}_vals']).split(',') if v.strip() not in ('', 'nan')]
        for v in vals:
            rows.append({'pair': label, 'group': g, 'ratio': v})
long_df = pd.DataFrame(rows)

n = len(group_names)
step = 0.8 / n
offsets = {g: (i - (n - 1) / 2) * step for i, g in enumerate(group_names)}

fig, ax = plt.subplots(figsize=(max(5, len(df) * 2.5), 5))
palette = sns.color_palette("Set2", n)
group_palette = dict(zip(group_names, palette))

sns.boxplot(data=long_df, x='pair', y='ratio', hue='group', ax=ax,
            palette=group_palette, showfliers=False, width=0.8, gap=0.12, legend=True)
sns.stripplot(data=long_df, x='pair', y='ratio', hue='group',
              dodge=True, color='black', alpha=0.7, size=4, ax=ax, legend=False)

y_min, y_max_plot = ax.get_ylim()
y_range = y_max_plot - y_min
bracket_step = y_range * 0.12

for x_idx, (_, row) in enumerate(df.iterrows()):
    for bracket_idx, (ga, gb) in enumerate(combinations(group_names, 2)):
        col = f'{ga}_{gb}_pvalue'
        if col not in df.columns:
            continue
        pval = row[col]
        x1 = x_idx + offsets[ga]
        x2 = x_idx + offsets[gb]
        y = y_max_plot + y_range * 0.05 + bracket_idx * bracket_step
        ax.plot([x1, x1, x2, x2], [y, y + y_range * 0.02, y + y_range * 0.02, y], 'k-', linewidth=0.8)
        label = f"p={pval:.2e}" if pval < 0.001 else f"p={pval:.3f}"
        ax.text((x1 + x2) / 2, y + y_range * 0.03, label, ha='center', fontsize=7,
                color='red' if pval < 0.05 else 'black')

n_pairs = len(list(combinations(group_names, 2)))
ax.set_ylim(y_min, y_max_plot + y_range * 0.05 + n_pairs * bracket_step + y_range * 0.05)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('')
ax.set_ylabel('Spatial ratio')
ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='right', fontsize=8)
ax.tick_params(axis='x', pad=5)
plt.tight_layout()
plt.savefig(args.output, dpi=330, transparent=False)
print(f"Saved plot to {args.output}")
