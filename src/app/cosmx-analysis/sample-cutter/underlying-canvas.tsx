

import React, { useEffect, useRef } from "react";

// HARD-CODING: FOR TESTING ONLY
import LarcBSampleDataCSV from "./larc-b-sample-data-csv";

// Efficient flat Float32Array: [x0, y0, x1, y1, ...]
const rawPoints = LarcBSampleDataCSV.split("\n").slice(1);
const positions = new Float32Array(rawPoints.length * 2);
rawPoints.forEach((line, i) => {
    const [, x, y] = line.split(",");
    positions[i * 2]     = parseFloat(x);
    positions[i * 2 + 1] = parseFloat(y);
});

const VERT = `
attribute vec2 a_pos;
uniform vec2 u_scale;  // 2 / (dataMax - dataMin)
uniform vec2 u_offset; // -1 - dataMin * scale
void main() {
    gl_Position = vec4(a_pos * u_scale + u_offset, 0, 1);
    gl_PointSize = 2.0;
}`;
const FRAG = `precision mediump float;
void main() { gl_FragColor = vec4(0.2, 0.5, 1.0, 0.7); }`;

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src); gl.compileShader(s);
    return s;
}

export function UnderlyingCanvas() {
    const canvasRef = useRef<HTMLCanvasElement>(null);

    useEffect(() => {
        const canvas = canvasRef.current!;
        const gl = canvas.getContext("webgl")!;

        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER, VERT));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, FRAG));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        // Upload all points once to GPU
        const buf = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, buf);
        gl.bufferData(gl.ARRAY_BUFFER, positions, gl.STATIC_DRAW);

        const a_pos = gl.getAttribLocation(prog, "a_pos");
        gl.enableVertexAttribArray(a_pos);
        gl.vertexAttribPointer(a_pos, 2, gl.FLOAT, false, 0, 0);

        // Compute data bounds for normalization to clip space [-1, 1]
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        for (let i = 0; i < positions.length; i += 2) {
            minX = Math.min(minX, positions[i]);     maxX = Math.max(maxX, positions[i]);
            minY = Math.min(minY, positions[i + 1]); maxY = Math.max(maxY, positions[i + 1]);
        }
        const sx = 2 / (maxX - minX), sy = 2 / (maxY - minY);
        gl.uniform2f(gl.getUniformLocation(prog, "u_scale"),  sx, sy);
        gl.uniform2f(gl.getUniformLocation(prog, "u_offset"), -1 - minX * sx, -1 - minY * sy);

        gl.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);
        // GPU clips points outside [-1,1] clip space automatically
        gl.drawArrays(gl.POINTS, 0, positions.length / 2);
    }, []);

    return <canvas ref={canvasRef} width={500} height={500} />;
}