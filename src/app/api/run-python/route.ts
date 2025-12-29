import { spawn } from 'child_process';

export async function GET(request: Request) {
    const { searchParams } = new URL(request.url);

    const pythonPath = process.env.PYTHON_PATH || "/mnt/sisplockers/ryanr/miniconda3/envs/scimap/bin/python";
    const scriptPath = searchParams.get('script');

    if (!scriptPath) {
        return Response.json({ error: 'script parameter is required' }, { status: 400 });
    }

    let args: string[] = [];
    searchParams.forEach((value, key) => {
        if (key !== 'script') {
            args.push(`--${key}`, value);
        }
    });
    args = ['-u', scriptPath, ...args];

    // The output file name if any
    let outputFilename: string | null = null;
    if (searchParams.get('output')) {
        outputFilename = searchParams.get('output') as string;
    }
    

    const encoder = new TextEncoder();
    const stream = new ReadableStream({
        start(controller) {
            let isClosed = false;

            const safeEnqueue = (data: Uint8Array) => {
                if (!isClosed) {
                    controller.enqueue(data);
                }
            };

            const safeClose = () => {
                if (!isClosed) {
                    isClosed = true;
                    controller.close();
                }
            };

            const childProcess = spawn(pythonPath, args, {
                cwd: 'workspace',
                env: { ...process.env, PYTHONUNBUFFERED: '1' }
            });

            // Send acknowledgment
            safeEnqueue(encoder.encode(`data: ${JSON.stringify({
                type: 'acknowledgment',
                status: 'processing',
                message: 'Script started'
            })}\n\n`));

            childProcess.stdout.on('data', (data) => {
                safeEnqueue(encoder.encode(`data: ${JSON.stringify({
                    type: 'stdout',
                    content: data.toString()
                })}\n\n`));
            });

            childProcess.stderr.on('data', (data) => {
                safeEnqueue(encoder.encode(`data: ${JSON.stringify({
                    type: 'stderr',
                    content: data.toString()
                })}\n\n`));
            });

            childProcess.on('close', (code) => {
                safeEnqueue(encoder.encode(`data: ${JSON.stringify({
                    type: 'complete',
                    exitCode: code,
                    outputFilename: outputFilename
                })}\n\n`));
                safeClose();
            });

            childProcess.on('error', (error) => {
                safeEnqueue(encoder.encode(`data: ${JSON.stringify({
                    type: 'error',
                    message: error.message
                })}\n\n`));
                safeClose();
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
    })

}