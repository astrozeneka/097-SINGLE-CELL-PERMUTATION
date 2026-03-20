"use client";
import { useEffect, useRef } from "react";

// Used by AbstractViewerController in viewer-2d (future shader injection phase).
export interface ShaderInjection {
    vsDeclarations: string;
    vsMainEnd: string;
}

// ---------------------------------------------------------------------------
// ColorEncoder<T>
// Fully describes how to color each point. The canvas compiles colorGlsl into
// the fragment shader and builds the interleaved buffer from the attributes.
//
// Built-in always available in colorGlsl (no need to declare in attributes):
//   a_selected  float  1.0 = selected, 0.0 = not selected  (from selection VBO)
//
// hue2rgb(h: float) -> vec3 is also always in scope.
// ---------------------------------------------------------------------------
export interface ColorEncoder<T> {
    attributes: {
        name: string;          // must follow "a_*" convention
        size: 1 | 2 | 3 | 4;
        feed: (d: T) => number | number[];
    }[];
    uniforms: {
        name: string;
        type: "float" | "vec2" | "vec3" | "vec4";
        value: number | number[];
    }[];
    // Body of the color function. Must contain `return vec4(...)`.
    // May reference any declared attribute by its a_* name or any uniform.
    colorGlsl: string;
}

export interface Transform { x: number; y: number; scale: number; }

interface UnderlyingCanvasParams<T> {
    data: T[];
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;
    zAccessor?: (d: T) => number; // reserved — projected away for 2D display
    colorEncoder: ColorEncoder<T>;
    transform: Transform;
    size: { w: number; h: number }; // driven by consumer ResizeObserver
    // Per-point float (1 = colored, 0 = greyed). null/undefined = all selected.
    // Uploaded to a dedicated VBO — position/attribute data is never re-uploaded.
    selectionMask?: Float32Array | null;
}

// ---------------------------------------------------------------------------
// Shader builder — compiles colorEncoder into vertex + fragment shaders.
// The vertex shader passes each ColorEncoder attribute as a varying so the
// fragment shader can reference them by their original a_* names via aliases.
// ---------------------------------------------------------------------------
function glslType(size: 1 | 2 | 3 | 4) { return size === 1 ? "float" : `vec${size}`; }
function varyingName(attrName: string)   { return "v" + attrName.slice(1); } // "a_cluster" → "v_cluster"

function buildVertShader(enc: ColorEncoder<any>): string {
    const attrDecls    = enc.attributes.map(a => `attribute ${glslType(a.size)} ${a.name};`).join("\n");
    const varyingDecls = enc.attributes.map(a => `varying ${glslType(a.size)} ${varyingName(a.name)};`).join("\n");
    const uniformDecls = enc.uniforms.map(u => `uniform ${u.type} ${u.name};`).join("\n");
    const assigns      = enc.attributes.map(a => `    ${varyingName(a.name)} = ${a.name};`).join("\n");
    return `
attribute vec2 a_pos;
${attrDecls}
attribute float a_selected;
uniform vec2 u_scale;
uniform vec2 u_offset;
uniform float u_userScale;
uniform vec2 u_userTranslate;
${varyingDecls}
varying float v_selected;
${uniformDecls}
void main() {
    vec2 clip = a_pos * u_scale + u_offset;
    clip = clip * u_userScale + u_userTranslate;
    gl_Position = vec4(clip, 0, 1);
    gl_PointSize = max(1.0, 2.0 * u_userScale);
${assigns}
    v_selected = a_selected;
}`;
}

function buildFragShader(enc: ColorEncoder<any>): string {
    const varyingDecls = enc.attributes.map(a => `varying ${glslType(a.size)} ${varyingName(a.name)};`).join("\n");
    const uniformDecls = enc.uniforms.map(u => `uniform ${u.type} ${u.name};`).join("\n");
    // Declare a_* aliases so colorGlsl can use attribute names directly (not v_*)
    const aliases = enc.attributes.map(a => `    ${glslType(a.size)} ${a.name} = ${varyingName(a.name)};`).join("\n");
    return `precision mediump float;
${varyingDecls}
varying float v_selected;
${uniformDecls}
vec3 hue2rgb(float h) {
    vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
    vec3 p = abs(fract(vec3(h) + k.xyz) * 6.0 - k.www);
    return clamp(p - k.xxx, 0.0, 1.0);
}
vec4 _color() {
${aliases}
    float a_selected = v_selected;
    ${enc.colorGlsl}
}
void main() { gl_FragColor = _color(); }`;
}

// ---------------------------------------------------------------------------

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src);
    gl.compileShader(s);
    return s;
}

type GlState = {
    gl: WebGLRenderingContext;
    count: number;
    scaleLoc: WebGLUniformLocation;
    offsetLoc: WebGLUniformLocation;
    userScaleLoc: WebGLUniformLocation;
    userTranslateLoc: WebGLUniformLocation;
    selBuffer: WebGLBuffer; // dedicated VBO updated independently on selection change
    bounds: { minX: number; maxX: number; minY: number; maxY: number };
};

export function UnderlyingCanvas<T>({ data, xAccessor, yAccessor, colorEncoder, transform, size, selectionMask }: UnderlyingCanvasParams<T>) {
    const canvasRef  = useRef<HTMLCanvasElement>(null);
    const glRef      = useRef<GlState | null>(null);
    // Ref so the init effect can read the current mask without it being a dependency.
    const selMaskRef = useRef(selectionMask);
    selMaskRef.current = selectionMask;

    // Full re-init when colorEncoder changes: recompile shaders, rebuild buffer, setup attributes.
    useEffect(() => {
        const gl = canvasRef.current!.getContext("webgl")!;
        gl.clearColor(0, 0, 0, 1);

        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER,   buildVertShader(colorEncoder)));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, buildFragShader(colorEncoder)));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        // Build interleaved buffer: [x, y, attr…,  x, y, attr…, …]
        const floatsPerPoint = 2 + colorEncoder.attributes.reduce((s, a) => s + a.size, 0);
        const buffer = new Float32Array(data.length * floatsPerPoint);
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (let i = 0; i < data.length; i++) {
            let off = i * floatsPerPoint;
            const x = xAccessor(data[i]), y = yAccessor(data[i]);
            buffer[off++] = x; buffer[off++] = y;
            for (const attr of colorEncoder.attributes) {
                const v = attr.feed(data[i]);
                if (typeof v === "number") buffer[off++] = v;
                else for (const n of v as number[]) buffer[off++] = n;
            }
            if (x < minX) minX = x; if (x > maxX) maxX = x;
            if (y < minY) minY = y; if (y > maxY) maxY = y;
        }

        // Main VBO (STATIC — never updated after upload)
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

        // Upload encoder uniforms
        for (const u of colorEncoder.uniforms) {
            const loc = gl.getUniformLocation(prog, u.name)!;
            if (typeof u.value === "number") gl.uniform1f(loc, u.value);
            else if (u.type === "vec2") gl.uniform2fv(loc, u.value as number[]);
            else if (u.type === "vec3") gl.uniform3fv(loc, u.value as number[]);
            else                        gl.uniform4fv(loc, u.value as number[]);
        }

        // Selection VBO (DYNAMIC — updated independently; a_selected always available in colorGlsl)
        const selBuffer = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, selBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, selMaskRef.current ?? new Float32Array(data.length).fill(1), gl.DYNAMIC_DRAW);
        const selLoc = gl.getAttribLocation(prog, "a_selected");
        gl.enableVertexAttribArray(selLoc);
        gl.vertexAttribPointer(selLoc, 1, gl.FLOAT, false, 0, 0);

        glRef.current = {
            gl, count: data.length,
            scaleLoc:         gl.getUniformLocation(prog, "u_scale")!,
            offsetLoc:        gl.getUniformLocation(prog, "u_offset")!,
            userScaleLoc:     gl.getUniformLocation(prog, "u_userScale")!,
            userTranslateLoc: gl.getUniformLocation(prog, "u_userTranslate")!,
            selBuffer,
            bounds: { minX, maxX, minY, maxY },
        };
    }, [colorEncoder]);

    // Upload selection mask to dedicated VBO only — main buffer untouched.
    // Declared before the draw effect so it always runs first in the same cycle.
    useEffect(() => {
        const state = glRef.current;
        if (!state) return;
        const { gl, selBuffer, count } = state;
        gl.bindBuffer(gl.ARRAY_BUFFER, selBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, selectionMask ?? new Float32Array(count).fill(1), gl.DYNAMIC_DRAW);
    }, [selectionMask]);

    // Redraw whenever size, transform, selection, or encoder changes.
    useEffect(() => {
        const state = glRef.current;
        if (!state || size.w === 0 || size.h === 0) return;
        const { gl, count, scaleLoc, offsetLoc, userScaleLoc, userTranslateLoc, bounds } = state;
        const { minX, maxX, minY, maxY } = bounds;
        const canvas = canvasRef.current!;

        if (canvas.width !== size.w || canvas.height !== size.h) {
            canvas.width  = size.w;
            canvas.height = size.h;
            gl.viewport(0, 0, size.w, size.h);
        }

        const dataW = maxX - minX, dataH = maxY - minY;
        const s  = (dataW / dataH > size.w / size.h) ? size.w / dataW : size.h / dataH;
        const sx = s * 2 / size.w, sy = s * 2 / size.h;
        gl.uniform2f(scaleLoc,  sx, sy);
        gl.uniform2f(offsetLoc, -(minX + dataW / 2) * sx, -(minY + dataH / 2) * sy);

        gl.uniform1f(userScaleLoc, transform.scale);
        gl.uniform2f(userTranslateLoc, 2 * transform.x / size.w, -2 * transform.y / size.h);

        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, count);
    }, [size, transform, selectionMask, colorEncoder]);

    return <canvas ref={canvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />;
}
