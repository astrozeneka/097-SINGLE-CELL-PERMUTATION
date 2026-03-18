import argparse
import pandas as pd
from itertools import combinations
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str, default='mannwhitney.csv')
args, unknown = parser.parse_known_args()

groups = {}
group_names = {}
current_key = None
i = 0
while i < len(unknown):
    arg = unknown[i]
    if arg.startswith('--groupname'):
        group_names[arg[2:].replace('name', '')] = unknown[i + 1] if i + 1 < len(unknown) else ''
        i += 2
    elif arg.startswith('--group'):
        current_key = arg[2:]
        groups.setdefault(current_key, [])
        i += 1
    elif current_key and not arg.startswith('--'):
        groups[current_key].append(arg)
        i += 1
    else:
        i += 1

group_keys = sorted(groups.keys())
label = lambda k: group_names.get(k, k)

def load_group(files):
    rows = {}
    for f in files:
        df = pd.read_csv(f)
        ratio_col = next(c for c in df.columns if c.startswith('d('))
        for _, row in df.iterrows():
            key = (row['target'], row['effector'])
            rows.setdefault(key, []).append(row[ratio_col])
    return rows

group_data = {k: load_group(v) for k, v in groups.items()}

all_pairs = set()
for data in group_data.values():
    all_pairs.update(data.keys())

records = []
for pair in sorted(all_pairs):
    group_vals = {k: group_data[k].get(pair, []) for k in group_keys}
    if any(len([v for v in group_vals[k] if pd.notna(v)]) < 2 for k in group_keys):
        continue
    record = {'target': pair[0], 'effector': pair[1]}
    for k in group_keys:
        vals = group_vals[k]
        record[f'{label(k)}_vals'] = ','.join([str(v) for v in vals if pd.notna(v)])
        record[f'{label(k)}_mean'] = sum([v for v in vals if pd.notna(v)]) / len([v for v in vals if pd.notna(v)])
    for ga, gb in combinations(group_keys, 2):
        stat, pval = mannwhitneyu([v for v in group_vals[ga] if pd.notna(v)], [v for v in group_vals[gb] if pd.notna(v)], alternative='two-sided')
        record[f'{label(ga)}_{label(gb)}_utest'] = stat
        record[f'{label(ga)}_{label(gb)}_pvalue'] = pval
    records.append(record)

pd.DataFrame(records).to_csv(args.output, index=False)
print(f"Saved {len(records)} rows to {args.output}")
