"use client";

import { useEffect, useRef } from "react";
import { Transform } from "./sample-cutter-client";

// Interleaved buffer: [x, y, cluster,  x, y, cluster, ...]
// stride = 3 floats = 12 bytes

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
    gl_PointSize = 0.5 * u_userScale;
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
    gl_FragColor = vec4(hue2rgb(mod(v_cluster, u_numClusters) / u_numClusters), 0.8);
}`;

function compileShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src); gl.compileShader(s);
    return s;
}

type GlState = {
    gl: WebGLRenderingContext;
    count: number;
    userScaleLoc: WebGLUniformLocation;
    userTranslateLoc: WebGLUniformLocation;
};

export function UnderlyingCanvas({ data, numClusters, transform }: { data: number[]; numClusters: number; transform: Transform }) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const glRef = useRef<GlState | null>(null);

    useEffect(() => {
        const canvas = canvasRef.current!;
        canvas.width = canvas.clientWidth;
        canvas.height = canvas.clientHeight;

        const gl = canvas.getContext("webgl")!;
        gl.viewport(0, 0, canvas.width, canvas.height);
        gl.clearColor(0, 0, 0, 1);
        const prog = gl.createProgram()!;
        gl.attachShader(prog, compileShader(gl, gl.VERTEX_SHADER, VERT));
        gl.attachShader(prog, compileShader(gl, gl.FRAGMENT_SHADER, FRAG));
        gl.linkProgram(prog);
        gl.useProgram(prog);

        const buf_data = new Float32Array(data);
        gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
        gl.bufferData(gl.ARRAY_BUFFER, buf_data, gl.STATIC_DRAW);

        const stride = 12;
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
        const midX = (minX + maxX) / 2;
        const midY = (minY + maxY) / 2;
        const s = Math.min(canvas.width, canvas.height) / (maxX - minX);
        const scaleX = s * 2 / canvas.width, scaleY = s * 2 / canvas.height;
        gl.uniform2f(gl.getUniformLocation(prog, "u_scale"), scaleX, scaleY);
        gl.uniform2f(gl.getUniformLocation(prog, "u_offset"), -midX * scaleX, -midY * scaleY);

        gl.uniform1f(gl.getUniformLocation(prog, "u_numClusters"), numClusters);

        glRef.current = {
            gl,
            count: buf_data.length / 3,
            userScaleLoc: gl.getUniformLocation(prog, "u_userScale")!,
            userTranslateLoc: gl.getUniformLocation(prog, "u_userTranslate")!,
        };
    }, []);

    useEffect(() => {
        const state = glRef.current;
        if (!state) return;
        const { gl, count, userScaleLoc, userTranslateLoc } = state;
        const { width, height } = canvasRef.current!;
        const s = transform.scale;
        gl.uniform1f(userScaleLoc, s);
        gl.uniform2f(userTranslateLoc,
            (s - 1) + 2 * transform.x / width,
            (1 - s) - 2 * transform.y / height,
        );
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, count);
    }, [transform]);

    return <canvas ref={canvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />;
}
