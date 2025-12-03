import { readFile } from 'fs/promises';
import { NextRequest } from 'next/server';
import path from 'path';

export async function GET(request: NextRequest) {
  try {
    const { searchParams } = new URL(request.url);
    const filename = searchParams.get('filename');

    if (!filename) {
      return Response.json({ error: 'filename parameter is required' }, { status: 400 });
    }

    // Check if it's a zip file (in workspace root) or csv file (in data folder)
    const isZipFile = filename.endsWith('.zip');
    const filePath = isZipFile
      ? path.join(process.cwd(), 'workspace', filename)
      : path.join(process.cwd(), 'workspace', 'data', filename);

    const fileBuffer = await readFile(filePath);

    const contentType = isZipFile ? 'application/zip' : 'text/csv';

    return new Response(new Uint8Array(fileBuffer), {
      headers: {
        'Content-Type': contentType,
        'Content-Disposition': `attachment; filename="${filename}"`,
      },
    });
  } catch (error: any) {
    return Response.json({ error: error.message }, { status: 500 });
  }
}