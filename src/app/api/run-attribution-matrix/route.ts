import { spawn } from 'child_process';

export async function GET(request: Request) {
  const { searchParams } = new URL(request.url);
  const filename = searchParams.get('filename');
  const markers = searchParams.get('markers');
  const phenotypes = searchParams.get('phenotypes');
  const phenotype = searchParams.get('phenotype') || 'Phenotype';

  if (!filename) {
    return Response.json({ error: 'filename parameter is required' }, { status: 400 });
  }

  if (!markers) {
    return Response.json({ error: 'markers parameter is required' }, { status: 400 });
  }

  if (!phenotypes) {
    return Response.json({ error: 'phenotypes parameter is required' }, { status: 400 });
  }

  const encoder = new TextEncoder();

  const stream = new ReadableStream({
    start(controller) {
      const pythonPath = "python";
      const scriptPath = "002_attribution_matrix.py";
      const outputFilename = `attribution_${filename.replace('.csv', '')}.png`;
      const args = [
        "--input", `data/${filename}`,
        "--markers", markers,
        "--phenotypes", phenotypes,
        "--phenotype", phenotype,
        "--output", `data/${outputFilename}`
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
          exitCode: code,
          outputFilename: outputFilename
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