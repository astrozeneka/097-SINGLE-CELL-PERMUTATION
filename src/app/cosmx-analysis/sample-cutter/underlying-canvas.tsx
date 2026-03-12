"use client";

import { useEffect, useRef } from "react";

// Interleaved buffer: [x, y, cluster,  x, y, cluster, ...]
// stride = 3 floats = 12 bytes

const VERT = `
attribute vec2 a_pos;
attribute float a_cluster;
uniform vec2 u_scale;
uniform vec2 u_offset;
varying float v_cluster;
void main() {
    gl_Position = vec4(a_pos * u_scale + u_offset, 0, 1);
    gl_PointSize = 2.0;
    v_cluster = a_cluster;
}`;

const FRAG = `precision mediump float;
varying float v_cluster;
vec3 hue2rgb(float h) {
    vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
    vec3 p = abs(fract(vec3(h) + k.xyz) * 6.0 - k.www);
    return clamp(p - k.xxx, 0.0, 1.0);
}
void main() {
    gl_FragColor = vec4(hue2rgb(mod(v_cluster, 20.0) / 20.0), 0.8);
}`;

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src); gl.compileShader(s);
    return s;
}

export function UnderlyingCanvas({ data }: { data: number[] }) {
    const canvasRef = useRef<HTMLCanvasElement>(null);

    useEffect(() => {
        const buf_data = new Float32Array(data); // [x, y, cluster, ...]
        const n = buf_data.length / 3;
        const canvas = canvasRef.current!;
        const gl = canvas.getContext("webgl")!;

        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER, VERT));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, FRAG));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, buf_data, gl.STATIC_DRAW);

        const stride = 12; // 3 floats * 4 bytes
        const a_pos = gl.getAttribLocation(prog, "a_pos");
        gl.enableVertexAttribArray(a_pos);
        gl.vertexAttribPointer(a_pos, 2, gl.FLOAT, false, stride, 0);

        const a_cluster = gl.getAttribLocation(prog, "a_cluster");
        gl.enableVertexAttribArray(a_cluster);
        gl.vertexAttribPointer(a_cluster, 1, gl.FLOAT, false, stride, 8);

        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (let i = 0; i < buf_data.length; i += 3) {
            minX = Math.min(minX, buf_data[i]);     maxX = Math.max(maxX, buf_data[i]);
            minY = Math.min(minY, buf_data[i + 1]); maxY = Math.max(maxY, buf_data[i + 1]);
        }
        const sx = 2 / (maxX - minX), sy = 2 / (maxY - minY);
        gl.uniform2f(gl.getUniformLocation(prog, "u_scale"),  sx, sy);
        gl.uniform2f(gl.getUniformLocation(prog, "u_offset"), -1 - minX * sx, -1 - minY * sy);

        gl.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, n);
    }, [data]);

    return <canvas ref={canvasRef} width={500} height={500} />;
}
