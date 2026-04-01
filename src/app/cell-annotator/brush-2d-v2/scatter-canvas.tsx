"use client";
import { useEffect, useRef } from "react";
import { ColorEncoder, Transform } from "../underlying-canvas";
import { Polygon } from "./polygon-manager-canvas";

interface ScatterCanvasProps<T> {
    data: T[];
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;
    colorEncoder: ColorEncoder<T>;
    transform: Transform;
    setTransform?: (t: Transform) => void;
    size: { w: number; h: number };
    // Per-point float key uploaded to GPU. 0 = no polygon.
    polygonMask?: Float32Array | null;
    // CPU-side lookup: float key → set of polygons owning that key.
    polygonMap?: Map<number, Set<Polygon>>;
    dotSize?: number;
    onReady?: () => void;
}

function fitTransform<T>(data: T[], xAccessor: (d: T) => number, yAccessor: (d: T) => number, w: number, h: number): Transform {
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

function glslType(size: 1 | 2 | 3 | 4) { return size === 1 ? "float" : `vec${size}`; }
function varyingName(attrName: string)   { return "v" + attrName.slice(1); }

function buildVertShader(enc: ColorEncoder<any>): string {
    const attrDecls    = enc.attributes.map(a => `attribute ${glslType(a.size)} ${a.name};`).join("\n");
    const varyingDecls = enc.attributes.map(a => `varying ${glslType(a.size)} ${varyingName(a.name)};`).join("\n");
    const uniformDecls = enc.uniforms.map(u => `uniform ${u.type} ${u.name};`).join("\n");
    const assigns      = enc.attributes.map(a => `    ${varyingName(a.name)} = ${a.name};`).join("\n");
    return `
attribute vec2 a_pos;
${attrDecls}
attribute float a_polygon;
uniform float u_pixelScale;
uniform vec2 u_pixelOffset;
uniform vec2 u_size;
uniform float u_dotSize;
uniform float u_baseScale;
${varyingDecls}
varying float v_polygon;
${uniformDecls}
void main() {
    float sx = a_pos.x * u_pixelScale + u_pixelOffset.x;
    float sy = -a_pos.y * u_pixelScale + u_pixelOffset.y;
    gl_Position = vec4(2.0 * sx / u_size.x - 1.0, 1.0 - 2.0 * sy / u_size.y, 0.0, 1.0);
    gl_PointSize = u_dotSize * (u_pixelScale / u_baseScale);
${assigns}
    v_polygon = a_polygon;
}`;
}

function buildFragShader(enc: ColorEncoder<any>): string {
    const varyingDecls = enc.attributes.map(a => `varying ${glslType(a.size)} ${varyingName(a.name)};`).join("\n");
    const uniformDecls = enc.uniforms.map(u => `uniform ${u.type} ${u.name};`).join("\n");
    const aliases = enc.attributes.map(a => `    ${glslType(a.size)} ${a.name} = ${varyingName(a.name)};`).join("\n");
    return `precision mediump float;
${varyingDecls}
varying float v_polygon;
${uniformDecls}
vec3 hue2rgb(float h) {
    vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
    vec3 p = abs(fract(vec3(h) + k.xyz) * 6.0 - k.www);
    return clamp(p - k.xxx, 0.0, 1.0);
}
vec4 _color() {
${aliases}
    float a_polygon = v_polygon;
    ${enc.colorGlsl}
}
void main() { gl_FragColor = _color(); }`;
}

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src);
    gl.compileShader(s);
    return s;
}

type GlState = {
    gl: WebGLRenderingContext;
    count: number;
    pixelScaleLoc: WebGLUniformLocation;
    pixelOffsetLoc: WebGLUniformLocation;
    sizeLoc: WebGLUniformLocation;
    dotSizeLoc: WebGLUniformLocation;
    baseScaleLoc: WebGLUniformLocation;
    polygonBuffer: WebGLBuffer;
};

export function ScatterCanvas<T>({ data, xAccessor, yAccessor, colorEncoder, transform, setTransform, size, polygonMask, polygonMap, dotSize = 2, onReady }: ScatterCanvasProps<T>) {
    const canvasRef      = useRef<HTMLCanvasElement>(null);
    const glRef          = useRef<GlState | null>(null);
    const polyMaskRef    = useRef(polygonMask);
    polyMaskRef.current  = polygonMask;
    const onReadyRef     = useRef(onReady);
    onReadyRef.current   = onReady;
    const readyCalledRef = useRef(false);
    const fittedDataRef     = useRef<T[] | null>(null);
    const baseScaleRef      = useRef(1);
    const setTransformRef   = useRef(setTransform);
    setTransformRef.current = setTransform;

    useEffect(() => {
        console.log("data", data);
        
        readyCalledRef.current = false;
        const gl = canvasRef.current!.getContext("webgl")!;
        gl.clearColor(0, 0, 0, 1);

        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER,   buildVertShader(colorEncoder)));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, buildFragShader(colorEncoder)));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        const floatsPerPoint = 2 + colorEncoder.attributes.reduce((s, a) => s + a.size, 0);
        const buffer = new Float32Array(data.length * floatsPerPoint);
        for (let i = 0; i < data.length; i++) {
            let off = i * floatsPerPoint;
            const x = xAccessor(data[i]), y = yAccessor(data[i]);
            buffer[off++] = x; buffer[off++] = y;
            for (const attr of colorEncoder.attributes) {
                const v = attr.feed(data[i]);
                if (typeof v === "number") buffer[off++] = v;
                else for (const n of v as number[]) buffer[off++] = n;
            }
        }

        const mainBuf = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, mainBuf);
        gl.bufferData(gl.ARRAY_BUFFER, buffer, gl.STATIC_DRAW);

        const stride = floatsPerPoint * 4;
        const posLoc = gl.getAttribLocation(prog, "a_pos");
        gl.enableVertexAttribArray(posLoc);
        gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, stride, 0);
        let attrByteOffset = 8;
        for (const attr of colorEncoder.attributes) {
            const loc = gl.getAttribLocation(prog, attr.name);
            gl.enableVertexAttribArray(loc);
            gl.vertexAttribPointer(loc, attr.size, gl.FLOAT, false, stride, attrByteOffset);
            attrByteOffset += attr.size * 4;
        }

        for (const u of colorEncoder.uniforms) {
            const loc = gl.getUniformLocation(prog, u.name)!;
            if (typeof u.value === "number") gl.uniform1f(loc, u.value);
            else if (u.type === "vec2") gl.uniform2fv(loc, u.value as number[]);
            else if (u.type === "vec3") gl.uniform3fv(loc, u.value as number[]);
            else                        gl.uniform4fv(loc, u.value as number[]);
        }
        const polygonBuffer = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, polygonBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, polyMaskRef.current ?? new Float32Array(data.length).fill(0), gl.DYNAMIC_DRAW);
        const polyLoc = gl.getAttribLocation(prog, "a_polygon");
        gl.enableVertexAttribArray(polyLoc);
        gl.vertexAttribPointer(polyLoc, 1, gl.FLOAT, false, 0, 0);

        glRef.current = {
            gl, count: data.length,
            pixelScaleLoc:  gl.getUniformLocation(prog, "u_pixelScale")!,
            pixelOffsetLoc: gl.getUniformLocation(prog, "u_pixelOffset")!,
            sizeLoc:        gl.getUniformLocation(prog, "u_size")!,
            dotSizeLoc:     gl.getUniformLocation(prog, "u_dotSize")!,
            baseScaleLoc:   gl.getUniformLocation(prog, "u_baseScale")!,
            polygonBuffer,
        };
    }, [data, colorEncoder]);

    useEffect(() => {
        const state = glRef.current;
        if (!state) return;
        const { gl, polygonBuffer, count } = state;
        gl.bindBuffer(gl.ARRAY_BUFFER, polygonBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, polygonMask ?? new Float32Array(count).fill(0), gl.DYNAMIC_DRAW);
    }, [polygonMask]);

    useEffect(() => {
        const state = glRef.current;
        if (!state || size.w === 0 || size.h === 0 || state.count === 0) return;
        if (setTransformRef.current && fittedDataRef.current !== data) {
            fittedDataRef.current = data;
            const fitted = fitTransform(data, xAccessor, yAccessor, size.w, size.h);
            baseScaleRef.current = fitted.scale;
            setTransformRef.current(fitted);
            return;
        }
        const { gl, count, pixelScaleLoc, pixelOffsetLoc, sizeLoc, dotSizeLoc, baseScaleLoc } = state;
        const canvas = canvasRef.current!;

        if (canvas.width !== size.w || canvas.height !== size.h) {
            canvas.width  = size.w;
            canvas.height = size.h;
            gl.viewport(0, 0, size.w, size.h);
        }

        gl.uniform1f(pixelScaleLoc, transform.scale);
        gl.uniform2f(pixelOffsetLoc, transform.x, transform.y);
        gl.uniform2f(sizeLoc, size.w, size.h);
        gl.uniform1f(dotSizeLoc, dotSize);
        gl.uniform1f(baseScaleLoc, baseScaleRef.current);

        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, count);
        if (!readyCalledRef.current) {
            readyCalledRef.current = true;
            onReadyRef.current?.();
        }
    }, [size, transform, polygonMask, colorEncoder, data, dotSize]);

    return <canvas ref={canvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />;
}
