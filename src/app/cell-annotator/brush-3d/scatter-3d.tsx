import { useEffect, useRef } from "react";


// --- sample data (25-point helix) ---
export const DATA = Array.from({ length: 25 }, (_, i) => {
    const t = (i / 24) * Math.PI * 4;
    return { x: Math.cos(t), y: (i / 24) * 2 - 1, z: Math.sin(t) };
});

// --- mat4 helpers (column-major, WebGL convention) ---
export function mat4mul(a: Float32Array, b: Float32Array): Float32Array {
    const out = new Float32Array(16);
    for (let j = 0; j < 4; j++)
        for (let i = 0; i < 4; i++)
            for (let k = 0; k < 4; k++)
                out[i + j * 4] += a[i + k * 4] * b[k + j * 4];
    return out;
}

export function rotateY(t: number): Float32Array {
    const c = Math.cos(t), s = Math.sin(t);
    return new Float32Array([c, 0, -s, 0,  0, 1, 0, 0,  s, 0, c, 0,  0, 0, 0, 1]);
}

export function rotateX(t: number): Float32Array {
    const c = Math.cos(t), s = Math.sin(t);
    return new Float32Array([1, 0, 0, 0,  0, c, s, 0,  0, -s, c, 0,  0, 0, 0, 1]);
}

export function perspectiveMat(fov: number, aspect: number, near: number, far: number): Float32Array {
    const f = 1 / Math.tan(fov / 2);
    const nf = 1 / (near - far);
    return new Float32Array([
        f / aspect, 0, 0,                    0,
        0,          f, 0,                    0,
        0,          0, (far + near) * nf,   -1,
        0,          0, 2 * far * near * nf,  0,
    ]);
}

// Camera at (0, 0, cameraZ) looking at origin along -Z
export function viewMat(cameraZ: number): Float32Array {
    return new Float32Array([1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, -cameraZ, 1]);
}


// --- WebGL shaders ---
const VS = `
attribute vec3 a_pos;
uniform mat4 u_matrix;
varying float v_depth;
void main() {
    vec4 p = u_matrix * vec4(a_pos, 1.0);
    gl_Position = p;
    gl_PointSize = 6.0;
    v_depth = p.w;
}`;

const FS = `
precision mediump float;
varying float v_depth;
void main() {
    float b = clamp(2.5 / v_depth, 0.2, 1.0);
    gl_FragColor = vec4(b, b * 0.5, 1.0, 1.0);
}`;

// --- Scatter3D (internal) ---
type GlState = { gl: WebGLRenderingContext; count: number; matrixLoc: WebGLUniformLocation };

export function Scatter3D({ matrix, size }: { matrix: Float32Array; size: { w: number; h: number } }) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const glRef     = useRef<GlState | null>(null);

    useEffect(() => {
        const gl = canvasRef.current!.getContext("webgl")!;
        gl.clearColor(0, 0, 0, 1);
        gl.enable(gl.DEPTH_TEST);

        const vs = gl.createShader(gl.VERTEX_SHADER)!;
        gl.shaderSource(vs, VS); gl.compileShader(vs);
        const fs = gl.createShader(gl.FRAGMENT_SHADER)!;
        gl.shaderSource(fs, FS); gl.compileShader(fs);
        const prog = gl.createProgram()!;
        gl.attachShader(prog, vs); gl.attachShader(prog, fs);
        gl.linkProgram(prog); gl.useProgram(prog);

        // normalize DATA to [-0.8, 0.8] using uniform scale
        let minX = Infinity, maxX = -Infinity;
        let minY = Infinity, maxY = -Infinity;
        let minZ = Infinity, maxZ = -Infinity;
        for (const d of DATA) {
            if (d.x < minX) minX = d.x; if (d.x > maxX) maxX = d.x;
            if (d.y < minY) minY = d.y; if (d.y > maxY) maxY = d.y;
            if (d.z < minZ) minZ = d.z; if (d.z > maxZ) maxZ = d.z;
        }
        const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ) / 1.6;
        const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2, cz = (minZ + maxZ) / 2;

        const buf = new Float32Array(DATA.length * 3);
        for (let i = 0; i < DATA.length; i++) {
            buf[i * 3]     = (DATA[i].x - cx) / span;
            buf[i * 3 + 1] = (DATA[i].y - cy) / span;
            buf[i * 3 + 2] = (DATA[i].z - cz) / span;
        }

        const vbo = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
        gl.bufferData(gl.ARRAY_BUFFER, buf, gl.STATIC_DRAW);
        const posLoc = gl.getAttribLocation(prog, "a_pos");
        gl.enableVertexAttribArray(posLoc);
        gl.vertexAttribPointer(posLoc, 3, gl.FLOAT, false, 0, 0);

        glRef.current = { gl, count: DATA.length, matrixLoc: gl.getUniformLocation(prog, "u_matrix")! };
    }, []);

    useEffect(() => {
        const state = glRef.current;
        if (!state || size.w === 0 || size.h === 0) return;
        const { gl, count, matrixLoc } = state;
        const canvas = canvasRef.current!;
        if (canvas.width !== size.w || canvas.height !== size.h) {
            canvas.width = size.w; canvas.height = size.h;
            gl.viewport(0, 0, size.w, size.h);
        }
        gl.uniformMatrix4fv(matrixLoc, false, matrix);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, count);
    }, [size, matrix]);

    return <canvas ref={canvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />;
}