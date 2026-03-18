import argparse
import pandas as pd
from itertools import combinations
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str, default='mannwhitney.csv')
args, unknown = parser.parse_known_args()

groups = {}
current_key = None
for arg in unknown:
    if arg.startswith('--group'):
        current_key = arg[2:]
        groups.setdefault(current_key, [])
    elif current_key:
        groups[current_key].append(arg)

group_keys = sorted(groups.keys())

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
    if any(len(group_vals[k]) < 2 for k in group_keys):
        continue
    record = {'target': pair[0], 'effector': pair[1]}
    for k in group_keys:
        vals = group_vals[k]
        record[f'{k}_vals'] = ','.join(str(v) for v in vals)
        record[f'{k}_mean'] = sum(vals) / len(vals)
    for ga, gb in combinations(group_keys, 2):
        stat, pval = mannwhitneyu(group_vals[ga], group_vals[gb], alternative='two-sided')
        record[f'{ga}_{gb}_utest'] = stat
        record[f'{ga}_{gb}_pvalue'] = pval
    records.append(record)

pd.DataFrame(records).to_csv(args.output, index=False)
print(f"Saved {len(records)} rows to {args.output}")
