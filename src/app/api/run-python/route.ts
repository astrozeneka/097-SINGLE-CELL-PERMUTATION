import { spawn } from 'child_process';

export async function GET(request: Request) {
    const { searchParams } = new URL(request.url);

    // const pythonPath = "/mnt/sisplockers/jantappapac/Ryan/conda/scimap/bin/python";
    const pythonPath = "D:/Ryan/126-NHOOD-WSA/.venv/Scripts/python.exe";
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
            // Send acknowledgment
            controller.enqueue(encoder.encode(`data: ${JSON.stringify({
                type: 'acknowledgment',
                status: 'processing',
                message: 'Script started'
            })}\n\n`));

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
    })

}