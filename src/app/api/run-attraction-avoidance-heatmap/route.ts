import { spawn } from 'child_process';
import { promisify } from 'util';
import { exec } from 'child_process';

const execAsync = promisify(exec);

export async function GET(request: Request) {
  const { searchParams } = new URL(request.url);
  const filenamesParam = searchParams.get('filenames');
  const displayDigitsParam = searchParams.get('displayDigits');

  if (!filenamesParam) {
    return Response.json({ error: 'filenames parameter is required' }, { status: 400 });
  }

  const encoder = new TextEncoder();

  const stream = new ReadableStream({
    start(controller) {
      const pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
      const scriptPath = "45b_attraction_avoidance_heatmap.py";

      // Parse filenames array and build selector pattern
      const filenames = JSON.parse(filenamesParam);

      // Generate unique output directory name using timestamp
      const timestamp = Date.now();
      const outputDir = `attraction_avoidance_heatmap_${timestamp}`;

      // Build selector pattern with proper escaping
      const selector = filenames.map((filename: string) => `data/${filename}`).join(' ');

      const args = [
        "--selector", selector,
        "--output-dir", outputDir
      ];

      if (displayDigitsParam === 'true') {
        args.push('--display-digits');
      }

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
