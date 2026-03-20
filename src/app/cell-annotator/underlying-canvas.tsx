"use client";
import { useEffect, useRef } from "react";

// Used by AbstractViewerController in viewer-2d (future shader injection phase).
export interface ShaderInjection {
    vsDeclarations: string;
    vsMainEnd: string;
}

// ---------------------------------------------------------------------------
// ColorEncoder<T>
// Describes how to color each point entirely from the consumer side.
// The canvas auto-generates shader declarations and builds the interleaved
// buffer — the consumer only writes the domain-specific parts.
//
// colorGlsl has access to every declared attribute and uniform by name.
// The hue2rgb() helper is always available in colorGlsl.
//
// Example (color by cluster):
//   const enc: ColorEncoder<Cell> = {
//       attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
//       uniforms:   [{ name: "u_n", type: "float", value: 14 }],
//       colorGlsl:  `return vec4(hue2rgb(a_cluster / u_n), 0.8);`,
//   };
// ---------------------------------------------------------------------------
export interface ColorEncoder<T> {
    // Extra per-point floats appended after x, y in the interleaved buffer.
    // Each attribute is auto-declared in the vertex shader and passed as a
    // varying to the fragment shader.
    attributes: {
        name: string;
        size: 1 | 2 | 3 | 4;
        feed: (d: T) => number | number[];
    }[];

    // Scalar/vector uniforms uploaded once before drawing.
    uniforms: {
        name: string;
        type: "float" | "vec2" | "vec3" | "vec4";
        value: number | number[];
    }[];

    // GLSL body of the color function — must contain a `return vec4(...)`.
    // May reference any declared attribute or uniform by name.
    // hue2rgb(h: float) -> vec3  is always available.
    colorGlsl: string;
}

interface UnderlyingCanvasParams<T> {
    data: T[];

    // Position accessors. The canvas auto-fits all points on init.
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;

    // Optional third dimension (e.g. t-SNE / UMAP z-axis).
    // Projected away for 2D display — reserved for future 3D support.
    zAccessor?: (d: T) => number;

    // Drives all coloring logic. Swap at runtime to recolor without
    // touching the canvas implementation (cluster → cell type → sample, etc.).
    colorEncoder: ColorEncoder<T>;

    // Controlled externally so sibling canvases stay in sync.
    transform: Transform;

    // Driven by a ResizeObserver in the consumer — triggers aspect ratio recompute.
    size: { w: number; h: number };
}

export interface Transform { x: number; y: number; scale: number; }

// Step 1: center + fit data → clip space  (u_scale encodes fit AND aspect ratio correction)
// Step 2: apply user zoom / pan
const VERT = `
attribute vec2 a_pos;
attribute float a_cluster;
uniform vec2 u_scale;
uniform vec2 u_offset;
uniform float u_userScale;
uniform vec2 u_userTranslate;
varying float v_cluster;
void main() {
    vec2 clip = a_pos * u_scale + u_offset;
    clip = clip * u_userScale + u_userTranslate;
    gl_Position = vec4(clip, 0, 1);
    gl_PointSize = max(1.0, 2.0 * u_userScale);
    v_cluster = a_cluster;
}`;

const FRAG = `precision mediump float;
varying float v_cluster;
uniform float u_numClusters;
vec3 hue2rgb(float h) {
    vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
    vec3 p = abs(fract(vec3(h) + k.xyz) * 6.0 - k.www);
    return clamp(p - k.xxx, 0.0, 1.0);
}
void main() {
    gl_FragColor = vec4(hue2rgb(v_cluster / u_numClusters), 0.8);
}`;

type GlState = {
    gl: WebGLRenderingContext;
    count: number;
    scaleLoc: WebGLUniformLocation;
    offsetLoc: WebGLUniformLocation;
    userScaleLoc: WebGLUniformLocation;
    userTranslateLoc: WebGLUniformLocation;
    bounds: { minX: number; maxX: number; minY: number; maxY: number };
};

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src);
    gl.compileShader(s);
    return s;
}

export function UnderlyingCanvas<T>({ data, xAccessor, yAccessor, colorEncoder, transform, size }: UnderlyingCanvasParams<T>) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const glRef = useRef<GlState | null>(null);

    // Init once: compile shaders, upload buffer, compute data bounds.
    // colorEncoder.attributes[0].feed used as cluster accessor until shader injection is implemented.
    useEffect(() => {
        const gl = canvasRef.current!.getContext("webgl")!;
        gl.clearColor(0, 0, 0, 1);

        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER, VERT));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, FRAG));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        const clusterFeed = colorEncoder.attributes[0].feed;
        let numClusters = 0;
        const raw: number[] = [];
        for (const d of data) {
            const c = clusterFeed(d) as number;
            raw.push(xAccessor(d), yAccessor(d), c);
            if (c > numClusters) numClusters = c;
        }
        const buffer = new Float32Array(raw);

        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, buffer, gl.STATIC_DRAW);

        const stride = 12; // 3 floats × 4 bytes
        const posLoc = gl.getAttribLocation(prog, "a_pos");
        gl.enableVertexAttribArray(posLoc);
        gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, stride, 0);

        const clusterLoc = gl.getAttribLocation(prog, "a_cluster");
        gl.enableVertexAttribArray(clusterLoc);
        gl.vertexAttribPointer(clusterLoc, 1, gl.FLOAT, false, stride, 8);

        gl.uniform1f(gl.getUniformLocation(prog, "u_numClusters")!, numClusters + 1);

        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (let i = 0; i < buffer.length; i += 3) {
            if (buffer[i]     < minX) minX = buffer[i];
            if (buffer[i]     > maxX) maxX = buffer[i];
            if (buffer[i + 1] < minY) minY = buffer[i + 1];
            if (buffer[i + 1] > maxY) maxY = buffer[i + 1];
        }

        glRef.current = {
            gl, count: buffer.length / 3,
            scaleLoc:         gl.getUniformLocation(prog, "u_scale")!,
            offsetLoc:        gl.getUniformLocation(prog, "u_offset")!,
            userScaleLoc:     gl.getUniformLocation(prog, "u_userScale")!,
            userTranslateLoc: gl.getUniformLocation(prog, "u_userTranslate")!,
            bounds: { minX, maxX, minY, maxY },
        };
    }, []);

    // Recompute aspect-correct fit + redraw on size or transform change.
    useEffect(() => {
        const state = glRef.current;
        if (!state || size.w === 0 || size.h === 0) return;
        const { gl, count, scaleLoc, offsetLoc, userScaleLoc, userTranslateLoc, bounds } = state;
        const { minX, maxX, minY, maxY } = bounds;
        const canvas = canvasRef.current!;

        // Only reset the drawing buffer when dimensions actually change.
        if (canvas.width !== size.w || canvas.height !== size.h) {
            canvas.width  = size.w;
            canvas.height = size.h;
            gl.viewport(0, 0, size.w, size.h);
        }

        // "Contain" fit: pick the axis that hits its limit first, so data never clips.
        const dataW = maxX - minX, dataH = maxY - minY;
        const s = (dataW / dataH > size.w / size.h) ? size.w / dataW : size.h / dataH;
        // ×2 converts pixel scale → clip range [-1, 1]. X and Y get separate factors
        // so the aspect ratio difference between data and canvas is absorbed here.
        const sx = s * 2 / size.w, sy = s * 2 / size.h;
        gl.uniform2f(scaleLoc,  sx, sy);
        gl.uniform2f(offsetLoc, -(minX + dataW / 2) * sx, -(minY + dataH / 2) * sy);

        // Pan is stored in screen pixels; converted to clip space here.
        // Y negated: screen Y increases downward, clip Y upward.
        gl.uniform1f(userScaleLoc, transform.scale);
        gl.uniform2f(userTranslateLoc,
             2 * transform.x / size.w,
            -2 * transform.y / size.h,
        );

        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, count);
    }, [size, transform]);

    return <canvas ref={canvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />;
}
