import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inputs")
parser.add_argument("--centroid_x_col")
parser.add_argument("--centroid_y_col")
parser.add_argument("--parent_area_col")
parser.add_argument("--cell_type_col")
args = parser.parse_args()

inputs = args.inputs.split(",") if args.inputs else []

print(f"inputs: {inputs}")
print(f"centroid_x_col: {args.centroid_x_col}")
print(f"centroid_y_col: {args.centroid_y_col}")
print(f"parent_area_col: {args.parent_area_col}")
print(f"cell_type_col: {args.cell_type_col}")
