import { exec } from 'child_process';
import { promisify } from 'util';

const execAsync = promisify(exec);

export async function GET() {
  try {
    let pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
    let scriptPath = "workspace/43_spatial_distance_permutation_test.py";
    let scriptParameter = " --input workspace/data/S19_12126B1_classified_cells_with_phenotypes.csv --n_permutations 10 --output workspace/data/permutation_results_S20_2317A9.csv --sample S20_2317A9";

    const command = `${pythonPath} ${scriptPath}${scriptParameter}`;
    const { stdout, stderr } = await execAsync(command);

    return Response.json({ output: stdout, error: stderr });
  } catch (error: any) {
    return Response.json({ error: error.message }, { status: 500 });
  }
}