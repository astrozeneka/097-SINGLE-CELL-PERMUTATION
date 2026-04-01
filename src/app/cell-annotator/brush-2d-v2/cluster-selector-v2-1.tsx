"use client";
import { useEffect, useState, useMemo } from "react";
import { ColorEncoder } from "../underlying-canvas";
import { CellDataV2_1 } from "./csv-upload-dialog-v2-1";

interface Props {
    data: CellDataV2_1[];
    colorEncoder: ColorEncoder<CellDataV2_1>;
    onMaskChange: (mask: Float32Array) => void;
}

function glslType(size: 1 | 2 | 3 | 4) { return size === 1 ? "float" : `vec${size}`; }

function buildSwatchShaders(enc: ColorEncoder<CellDataV2_1>) {
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

function renderSwatchToDataUrl(clusterIdx: number, enc: ColorEncoder<CellDataV2_1>): string {
    const canvas = document.createElement("canvas");
    canvas.width = 24; canvas.height = 24;
    const gl = canvas.getContext("webgl")!;
    const { vert, frag } = buildSwatchShaders(enc);

    const compile = (type: number, src: string) => {
        const s = gl.createShader(type)!;
        gl.shaderSource(s, src); gl.compileShader(s); return s;
    };
    const prog = gl.createProgram()!;
    gl.attachShader(prog, compile(gl.VERTEX_SHADER, vert));
    gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, frag));
    gl.linkProgram(prog); gl.useProgram(prog);

    const quad = new Float32Array([-1,-1, 1,-1, -1,1, -1,1, 1,-1, 1,1]);
    const posBuf = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.bufferData(gl.ARRAY_BUFFER, quad, gl.STATIC_DRAW);
    const posLoc = gl.getAttribLocation(prog, "a_pos");
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    const dummy = { clusterIdx } as unknown as CellDataV2_1;
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

export default function ClusterSelectorV2_1({ data, colorEncoder, onMaskChange }: Props) {
    const clusters = useMemo(() => {
        const seen = new Map<string, number>();
        for (const d of data) if (!seen.has(d.cluster)) seen.set(d.cluster, d.clusterIdx);
        return Array.from(seen.entries()).map(([name, idx]) => ({ name, idx })).sort((a, b) => a.idx - b.idx);
    }, [data]);

    const [swatches, setSwatches] = useState<Map<string, string>>(new Map());
    const [hidden, setHidden] = useState<Set<string>>(new Set());

    useEffect(() => {
        const map = new Map<string, string>();
        for (const { name, idx } of clusters) map.set(name, renderSwatchToDataUrl(idx, colorEncoder));
        setSwatches(map);
    }, [clusters, colorEncoder]);

    useEffect(() => {
        const mask = new Float32Array(data.length);
        for (let i = 0; i < data.length; i++) mask[i] = hidden.has(data[i].cluster) ? 0.0 : 1.0;
        onMaskChange(mask);
    }, [hidden, data]);

    const toggle = (name: string) =>
        setHidden(prev => { const s = new Set(prev); s.has(name) ? s.delete(name) : s.add(name); return s; });

    const selectAllState = hidden.size === 0 ? "checked" : hidden.size === clusters.length ? "unchecked" : "partial";
    const toggleAll = () => setHidden(selectAllState === "checked" ? new Set(clusters.map(c => c.name)) : new Set());

    return (
        <div className="thin-scrollbar" style={{
            position: "absolute", top: 8, right: 8, zIndex: 20,
            background: "rgba(0,0,0,0.55)", backdropFilter: "blur(4px)",
            padding: "6px 8px", display: "flex", flexDirection: "column", gap: 2,
            maxHeight: "80%", overflowY: "auto", fontSize: 11, fontFamily: "monospace",
        }}>
            <div onClick={toggleAll}
                style={{ display: "flex", alignItems: "center", gap: 6, cursor: "pointer", paddingBottom: 4, marginBottom: 2, borderBottom: "1px solid rgba(255,255,255,0.15)" }}>
                <div style={{
                    width: 14, height: 14, flexShrink: 0, border: "1.5px solid #aaa", borderRadius: 2,
                    display: "flex", alignItems: "center", justifyContent: "center",
                    background: selectAllState === "checked" ? "#aaa" : "transparent",
                }}>
                    {selectAllState === "partial" && <div style={{ width: 8, height: 2, background: "#aaa", borderRadius: 1 }} />}
                    {selectAllState === "checked" && <svg width="10" height="8" viewBox="0 0 10 8"><polyline points="1,4 4,7 9,1" fill="none" stroke="#222" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" /></svg>}
                </div>
                <span style={{ color: "#e2e8f0" }}>All</span>
            </div>
            {clusters.map(({ name, idx }) => {
                const isHidden = hidden.has(name);
                return (
                    <div key={name} onClick={() => toggle(name)}
                        style={{ display: "flex", alignItems: "center", gap: 6, cursor: "pointer", opacity: isHidden ? 0.35 : 1 }}>
                        {swatches.has(name)
                            ? <img src={swatches.get(name)} width={16} height={16} style={{ display: "block", flexShrink: 0 }} />
                            : <div style={{ width: 16, height: 16, background: "#333", flexShrink: 0 }} />
                        }
                        <span style={{ color: "#e2e8f0" }}>{name}</span>
                    </div>
                );
            })}
        </div>
    );
}
