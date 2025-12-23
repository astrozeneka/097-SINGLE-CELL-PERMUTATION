import { writeFile, mkdir } from 'fs/promises';
import { NextRequest } from 'next/server';
import path from 'path';
import crypto from 'crypto';

export async function POST(request: NextRequest) {
  try {
    const formData = await request.formData();
    const file = formData.get('file') as File;

    if (!file) {
      return Response.json({ error: 'No file provided' }, { status: 400 });
    }

    const bytes = await file.arrayBuffer();
    const buffer = Buffer.from(bytes);

    const uploadDir = path.join(process.cwd(), 'workspace', 'data');
    await mkdir(uploadDir, { recursive: true });

    const ext = path.extname(file.name);
    const nameWithoutExt = path.basename(file.name, ext);
    const randomValue = crypto.randomBytes(4).toString('hex');
    const newFilename = `${nameWithoutExt}_${randomValue}${ext}`;

    const filePath = path.join(uploadDir, newFilename);
    await writeFile(filePath, buffer);

    const relativePath = path.join('data', newFilename);

    return Response.json({
      success: true,
      filename: newFilename,
      path: relativePath
    });
  } catch (error: any) {
    return Response.json({ error: error.message }, { status: 500 });
  }
}