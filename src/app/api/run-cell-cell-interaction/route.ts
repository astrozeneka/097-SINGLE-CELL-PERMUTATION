import { spawn } from 'child_process';
import { promisify } from 'util';
import { exec } from 'child_process';

const execAsync = promisify(exec);

export async function GET(request: Request) {
  const { searchParams } = new URL(request.url);
  const fileTuplesParam = searchParams.get('file_tuples');
  const outputDir = searchParams.get('output_dir') || 'attraction_repulsion_results';

  if (!fileTuplesParam) {
    return Response.json({ error: 'file_tuples parameter is required' }, { status: 400 });
  }

  const encoder = new TextEncoder();

  const stream = new ReadableStream({
    start(controller) {
      const pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
      const scriptPath = "44_analyze_permutation_test.py";

      // Parse file tuples and convert to the format expected by the Python script
      const fileTuples = JSON.parse(fileTuplesParam);
      // Convert from [slug, observed_filename, permuted_filename] to [slug, observed_path, permuted_path]
      const fileTuplesWithPaths = fileTuples.map(([slug, observedFile, permutedFile]: [string, string, string]) => [
        slug,
        `data/${observedFile}`,
        `data/${permutedFile}`
      ]);

      const args = [
        "--file-tuple", JSON.stringify(fileTuplesWithPaths),
        "--output-dir", outputDir
      ];

      const childProcess = spawn(pythonPath, ['-u', scriptPath, ...args], {
        cwd: 'workspace',
        env: { ...process.env, PYTHONUNBUFFERED: '1' }
      });

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

      childProcess.on('close', async (code) => {
        if (code === 0) {
          // Zip the output directory
          try {
            controller.enqueue(encoder.encode(`data: ${JSON.stringify({
              type: 'stdout',
              content: 'Creating zip file...'
            })}\n\n`));

            const zipFilename = `${outputDir}.zip`;
            const zipPath = `workspace/${zipFilename}`;

            await execAsync(`cd workspace && zip -r ${zipFilename} ${outputDir}`);

            controller.enqueue(encoder.encode(`data: ${JSON.stringify({
              type: 'complete',
              exitCode: code,
              zipFilename: zipFilename
            })}\n\n`));
          } catch (zipError: any) {
            controller.enqueue(encoder.encode(`data: ${JSON.stringify({
              type: 'error',
              message: `Failed to create zip: ${zipError.message}`
            })}\n\n`));
          }
        } else {
          controller.enqueue(encoder.encode(`data: ${JSON.stringify({
            type: 'complete',
            exitCode: code
          })}\n\n`));
        }
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
