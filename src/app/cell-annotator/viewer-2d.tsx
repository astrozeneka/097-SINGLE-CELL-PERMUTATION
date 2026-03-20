"use client";
import { useEffect } from "react";
import sample_data_csv from "./brush-2d/sample_data_csv";
import { ShaderInjection, UnderlyingCanvas } from "./underlying-canvas";
import { OverlyingCanvas } from "./overlying-canvas";

interface _Cell {
    id: string;
    x: number;
    y: number;
    cluster: number;
}

interface Viewer2dParams {
}

abstract class AbstractViewerController<T> {
    abstract readonly injection: ShaderInjection;
    abstract buildBuffer(data: T[]): Float32Array;  // interleaved: [x, y, ...extras]

    abstract uploadUniforms(gl: WebGLRenderingContext, prog: WebGLProgram, data: T[]): void;
}

interface ColorEncoder {
    // Extra per-point floats this encoder needs (appended after x,y in buffer)
    attributes: { name: string; size: 1 | 2 | 3 | 4; feed: (d: any) => number | number[] }[];
    // Uniforms this encoder needs uploaded
    uniforms: { name: string; type: "float" | "vec2" | "vec3" | "vec4" }[];
    // GLSL function body: given the declared attributes/uniforms, return vec4
    // Only this expression is injected — declarations are auto-generated
    colorGlsl: string;
}

const byColorEncoder: ColorEncoder = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms:   [{ name: "u_numClusters", type: "float" }],
    colorGlsl:  `return vec4(hue2rgb(a_cluster / u_numClusters), 0.8);`,
}


export default function Viewer2d(params: Viewer2dParams) {
    useEffect(() => {
        // DRAFT CODE FOR TESTING
        const data: _Cell[] = sample_data_csv.split("\n").slice(1).map(line => {
            const [id, x, y, cluster] = line.split(",");
            return {
                id,
                x: parseFloat(x),
                y: parseFloat(y),
                z: 0,
                cluster: parseInt(cluster)
            };
        });
    }, []);

    return (
        <div>
            <OverlyingCanvas></OverlyingCanvas>
            <UnderlyingCanvas></UnderlyingCanvas>
        </div>
    );
}