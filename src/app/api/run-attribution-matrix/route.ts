/*
--input
"used_for_attribution_matrix/SAB_combined_zscore_838r.tif - 3 roi neighbourhood measurement_updated.csv"
--markers
CD68,CD163,CD206,CD11c,PDL1,CD3,CD4,Foxp3,CD8,PD1,CD20,CD31,Ki67
--phenotypes
"TAM,DC,T,Helper T,T PD1,Tregs,Endothelial cells,Ki67 proliferating cells"
--phenotype
Phenotype
--output
attribution.png
*/

export async function GET() {
  return Response.json({ message: "Hello World!" });
}