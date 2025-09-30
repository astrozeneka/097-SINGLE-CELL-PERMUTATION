import { exec } from 'child_process';
import { promisify } from 'util';

const execAsync = promisify(exec);

export async function GET() {
  try {
    const { stdout } = await execAsync('squeue');
    return Response.json({ output: stdout });
    let pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
    let scriptPath = "workspace/43_spatial_distance_permutation_test.py";
    let scriptParameter = " --selector workspace/data/cell_data_with_regions_and_parent_areas_S20_2317A9.csv --n_permutations 100"
    // Python path
    // /mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python
    /*
     {
    "output": "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)\n"
    */
  } catch (error: any) {
    return Response.json({ error: error.message }, { status: 500 });
  }
}