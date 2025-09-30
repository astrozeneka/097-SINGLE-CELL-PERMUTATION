import { exec } from 'child_process';
import { promisify } from 'util';

const execAsync = promisify(exec);

export async function GET() {
  try {
    const { stdout } = await execAsync('squeue');
    return Response.json({ output: stdout });
    /*
     {
    "output": "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)\n"
    */
  } catch (error: any) {
    return Response.json({ error: error.message }, { status: 500 });
  }
}