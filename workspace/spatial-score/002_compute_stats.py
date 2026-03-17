import argparse
import os
import shutil
import zipfile
import pandas as pd

parser = argparse.ArgumentParser(description="Compute spatial score statistics")
parser.add_argument('--reference-cell', type=str, required=True)
parser.add_argument('--target-cells', type=str, action='append', required=True)
parser.add_argument('--effector-cells', type=str, action='append', required=True)
parser.add_argument('--output', type=str, default='data/spatial_score_stats.zip')
args, unknown = parser.parse_known_args()

input_files = [v for k, v in zip(unknown[::2], unknown[1::2]) if k.startswith('--input')]

if __name__ == '__main__':
    dfs = []
    for f in input_files:
        print(f"Reading {f}...")
        dfs.append(pd.read_csv(f))
    df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(df)} cells total")

    tmp_dir = "spatial_score_raw"
    os.makedirs(tmp_dir, exist_ok=True)

    for image_id, image_df in df.groupby('imageid'):
        print(f"Processing imageid={image_id} ({len(image_df)} cells)...")
        ratio_rows = []
        for target in args.target_cells:
            sample_df = image_df[image_df['Phenotype'] == target].copy()
            if sample_df.empty:
                continue
            if args.reference_cell not in sample_df.columns:
                print(f"  WARNING: reference column '{args.reference_cell}' not found, skipping")
                continue
            for effector in args.effector_cells:
                if effector not in sample_df.columns:
                    print(f"  WARNING: effector column '{effector}' not found, skipping")
                    continue
                nearest_reference = sample_df[args.reference_cell]
                nearest_effector = sample_df[effector]
                ratio = nearest_reference.div(nearest_effector).where(nearest_effector > 0)
                ratio_rows.append({
                    'target': target,
                    'effector': effector,
                    f'd(target, {args.reference_cell}) / d(target, effector)': ratio.mean(),
                    'count': ratio.count(),
                    'imageid': image_id,
                })

        pd.DataFrame(ratio_rows).to_csv(f"{tmp_dir}/{image_id}.csv", index=False)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with zipfile.ZipFile(args.output, 'w') as z:
        for f in os.listdir(tmp_dir):
            z.write(f"{tmp_dir}/{f}", f)
    shutil.rmtree(tmp_dir)
    print(f"Output saved to {args.output}")
    print("Done")
