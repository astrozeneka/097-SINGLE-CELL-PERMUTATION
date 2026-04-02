import { Transform } from "../underlying-canvas";

export function glslType(size: 1 | 2 | 3 | 4) { return size === 1 ? "float" : `vec${size}`; }

export function varyingName(attrName: string)   { return "v" + attrName.slice(1); }

export function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src);
    gl.compileShader(s);
    return s;
}

export function fitTransform<T>(data: T[], xAccessor: (d: T) => number, yAccessor: (d: T) => number, w: number, h: number): Transform {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const d of data) {
        const x = xAccessor(d), y = yAccessor(d);
        if (x < minX) minX = x; if (x > maxX) maxX = x;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
    }
    const dataW = maxX - minX, dataH = maxY - minY;
    const scale = dataW / dataH > w / h ? w / dataW : h / dataH;
    const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
    return { x: w / 2 - cx * scale, y: h / 2 + cy * scale, scale };
}