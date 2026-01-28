import { spawn } from 'child_process';

export async function GET(request: Request) {
  const { searchParams } = new URL(request.url);
  const filename = searchParams.get('filename');
  const slug = searchParams.get('slug');
  const phenotype = searchParams.get('phenotype') || 'Phenotype';
  const xCoord = searchParams.get('x_coord') || 'Centroid X';
  const yCoord = searchParams.get('y_coord') || 'Centroid Y';

  if (!filename) {
    return Response.json({ error: 'filename parameter is required' }, { status: 400 });
  }

  if (!slug) {
    return Response.json({ error: 'slug parameter is required' }, { status: 400 });
  }

  const encoder = new TextEncoder();

  const stream = new ReadableStream({
    start(controller) {
      const pythonPath = "python";
      const scriptPath = "62_spatial_distance_real_tissue.py";
      const outputFilename = `spatial_distance_${filename}`;
      const args = [
        "--input", `data/${filename}`,
        "--output", `data/${outputFilename}`,
        "--slug", slug,
        "--phenotype", phenotype,
        "--x-coord", xCoord,
        "--y-coord", yCoord
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
