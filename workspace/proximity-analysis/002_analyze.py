import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--group1_inputs")
parser.add_argument("--group2_inputs")
parser.add_argument("--output")
args = parser.parse_args()

# sample code, print the inputs
print("Group 1 inputs:", args.group1_inputs)
print("Group 2 inputs:", args.group2_inputs)

# Create dummy output file
df = pd.DataFrame({
    "Phenotype A": ["B_Cells", "T_Cells", "Macrophages"],
    "Phenotype B": ["B_Cells", "T_Cells", "Macrophages"],
})
df.to_csv(args.output, index=False)