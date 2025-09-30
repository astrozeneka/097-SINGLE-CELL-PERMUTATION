import { spawn } from 'child_process';

export async function GET() {
  const encoder = new TextEncoder();

  const stream = new ReadableStream({
    start(controller) {
      const pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
      const scriptPath = "43_spatial_distance_permutation_test.py";
      const args = [
        "--input", "data/S19_12126B1_classified_cells_with_phenotypes.csv",
        "--n_permutations", "10",
        "--output", "data/permutation_results_S20_2317A9.csv",
        "--sample", "S20_2317A9"
      ];

      const childProcess = spawn(pythonPath, [scriptPath, ...args], { cwd: 'workspace' });

      // Send acknowledgment
      controller.enqueue(encoder.encode(`data: ${JSON.stringify({
        type: 'acknowledgment',
        status: 'processing',
        message: 'Script started'
      })}\n\n`));

      childProcess.stdout.on('data', (data) => {
        controller.enqueue(encoder.encode(`data: ${JSON.stringify({
          type: 'stdout',
          content: data.toString()
        })}\n\n`));
      });

      childProcess.stderr.on('data', (data) => {
        controller.enqueue(encoder.encode(`data: ${JSON.stringify({
          type: 'stderr',
          content: data.toString()
        })}\n\n`));
      });

      childProcess.on('close', (code) => {
        controller.enqueue(encoder.encode(`data: ${JSON.stringify({
          type: 'complete',
          exitCode: code
        })}\n\n`));
        controller.close();
      });

      childProcess.on('error', (error) => {
        controller.enqueue(encoder.encode(`data: ${JSON.stringify({
          type: 'error',
          message: error.message
        })}\n\n`));
        controller.close();
      });
    }
  });

  return new Response(stream, {
    headers: {
      'Cache-Control': 'no-cache',
      'Content-Type': 'text/event-stream',
      'Connection': 'keep-alive',
      'X-Accel-Buffering': 'no',
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Headers': 'Cache-Control'
    }
  });
}