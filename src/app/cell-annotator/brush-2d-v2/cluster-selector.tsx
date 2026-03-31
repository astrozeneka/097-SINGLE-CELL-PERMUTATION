"use client";
import { useEffect, useState, useMemo } from "react";
import { ColorEncoder } from "../underlying-canvas";
import { CellDataV2 } from "./csv-upload-dialog-v2";

interface Props {
    data: CellDataV2[];
    colorEncoder: ColorEncoder<CellDataV2>;
}

// ─── Swatch WebGL box ────────────────────────────────────────────────────────
function glslType(size: 1 | 2 | 3 | 4) { return size === 1 ? "float" : `vec${size}`; }

function buildSwatchShaders(enc: ColorEncoder<CellDataV2>) {
    const attrDecls    = enc.attributes.map(a => `attribute ${glslType(a.size)} ${a.name};`).join("\n");
    const varyingDecls = enc.attributes.map(a => `varying ${glslType(a.size)} v_${a.name.slice(2)};`).join("\n");
    const assigns      = enc.attributes.map(a => `    v_${a.name.slice(2)} = ${a.name};`).join("\n");
    const aliases      = enc.attributes.map(a => `    ${glslType(a.size)} ${a.name} = v_${a.name.slice(2)};`).join("\n");

    const vert = `
attribute vec2 a_pos;
${attrDecls}
${varyingDecls}
void main() {
    gl_Position = vec4(a_pos, 0.0, 1.0);
${assigns}
}`;

    const frag = `
precision mediump float;
${varyingDecls}
vec3 hue2rgb(float h) {
    vec4 k = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
    vec3 p = abs(fract(vec3(h) + k.xyz) * 6.0 - k.www);
    return clamp(p - k.xxx, 0.0, 1.0);
}
vec4 _color() {
${aliases}
    float a_polygon = 1.0;
    ${enc.colorGlsl}
}
void main() { gl_FragColor = _color(); }`;

    return { vert, frag };
}

function renderSwatchToDataUrl(clusterId: number, enc: ColorEncoder<CellDataV2>): string {
    const canvas = document.createElement("canvas");
    canvas.width = 24;
    canvas.height = 24;
    const gl = canvas.getContext("webgl")!;
    const { vert, frag } = buildSwatchShaders(enc);

    const compile = (type: number, src: string) => {
        const s = gl.createShader(type)!;
        gl.shaderSource(s, src);
        gl.compileShader(s);
        return s;
    };
    const prog = gl.createProgram()!;
    gl.attachShader(prog, compile(gl.VERTEX_SHADER, vert));
    gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, frag));
    gl.linkProgram(prog);
    gl.useProgram(prog);

    const quad = new Float32Array([-1, -1,  1, -1,  -1, 1,  -1, 1,  1, -1,  1, 1]);
    const posBuf = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.bufferData(gl.ARRAY_BUFFER, quad, gl.STATIC_DRAW);
    const posLoc = gl.getAttribLocation(prog, "a_pos");
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    const dummy = { cluster: clusterId } as CellDataV2;
    for (const attr of enc.attributes) {
        const val = attr.feed(dummy);
        const flat = typeof val === "number" ? [val] : (val as number[]);
        const arr = new Float32Array(6 * attr.size);
        for (let i = 0; i < 6; i++) flat.forEach((v, j) => { arr[i * attr.size + j] = v; });
        const buf = gl.createBuffer()!;
        gl.bindBuffer(gl.ARRAY_BUFFER, buf);
        gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW);
        const loc = gl.getAttribLocation(prog, attr.name);
        gl.enableVertexAttribArray(loc);
        gl.vertexAttribPointer(loc, attr.size, gl.FLOAT, false, 0, 0);
    }

    gl.drawArrays(gl.TRIANGLES, 0, 6);
    const dataUrl = canvas.toDataURL();
    gl.getExtension("WEBGL_lose_context")?.loseContext();
    return dataUrl;
}
// ─────────────────────────────────────────────────────────────────────────────

export default function ClusterSelector({ data, colorEncoder }: Props) {
    const clusters = useMemo(
        () => Array.from(new Set(data.map(d => d.cluster))).sort((a, b) => a - b),
        [data]
    );

    const [swatches, setSwatches] = useState<Map<number, string>>(new Map());
    const [hidden, setHidden] = useState<Set<number>>(new Set());

    useEffect(() => {
        const map = new Map<number, string>();
        for (const c of clusters) map.set(c, renderSwatchToDataUrl(c, colorEncoder));
        setSwatches(map);
    }, [clusters, colorEncoder]);

    const toggle = (c: number) =>
        setHidden(prev => { const s = new Set(prev); s.has(c) ? s.delete(c) : s.add(c); return s; });

    return (
        <div style={{
            position: "absolute", top: 8, right: 8, zIndex: 20,
            background: "rgba(0,0,0,0.55)", backdropFilter: "blur(4px)",
            padding: "6px 8px", display: "flex", flexDirection: "column", gap: 2,
            maxHeight: "80%", overflowY: "auto", fontSize: 11, fontFamily: "monospace",
        }}>
            {clusters.map(c => {
                const isHidden = hidden.has(c);
                return (
                    <div
                        key={c}
                        onClick={() => toggle(c)}
                        style={{ display: "flex", alignItems: "center", gap: 6, cursor: "pointer", opacity: isHidden ? 0.35 : 1 }}
                    >
                        {swatches.has(c)
                            ? <img src={swatches.get(c)} width={16} height={16} style={{ display: "block", flexShrink: 0 }} />
                            : <div style={{ width: 16, height: 16, background: "#333", flexShrink: 0 }} />
                        }
                        <span style={{ color: "#e2e8f0" }}>{c}</span>
                    </div>
                );
            })}
        </div>
    );
}
