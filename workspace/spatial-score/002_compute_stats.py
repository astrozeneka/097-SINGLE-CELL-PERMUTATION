
import argparse

parser = argparse.ArgumentParser(description="Compute spatial score statistics")
parser.add_argument('--reference-cell', type=str, required=True, help="Phenotype to use as reference cell type (generally Tumor cells)")
parser.add_argument('--target-cells', type=str, nargs='+', required=True, help="Phenotypes to compute spatial scores for (generally immune cell types)")
parser.add_argument('--effector-cells', type=str, nargs='+', required=True, help="Phenotypes to compute spatial scores for (generally immune cell types)")
parser.add_argument('--output', type=str, default='data/spatial_score_stats.csv')

# Parse known arguments
args, unknown = parser.parse_known_args()

# Collect all --inputN arguments
input_files = [v for k, v in zip(unknown[::2], unknown[1::2]) if k.startswith('--input')]

if __name__ == '__main__':
    # Print all arguments as a dictionary
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")
    
    # TODO Later