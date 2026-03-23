"use client";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import CsvUploadDialog, { ParsedCsvData } from "./csv-upload-dialog";
import { ColorEncoder, ShaderInjection, Transform, UnderlyingCanvas } from "./underlying-canvas";
import { CanvasMode, OverlyingCanvas } from "./overlying-canvas";

type _Cell = ParsedCsvData["data"][number];

interface Viewer2dParams {
}

abstract class AbstractViewerController<T> {
    abstract readonly injection: ShaderInjection;
    abstract buildBuffer(data: T[]): Float32Array;  // interleaved: [x, y, ...extras]

    abstract uploadUniforms(gl: WebGLRenderingContext, prog: WebGLProgram, data: T[]): void;
}

// Encoders are defined inside the component (see useMemo below) so numClusters
// can be included as a uniform value once parsed from the data.

// Screen pixel → data coordinate. Mirrors the vertex shader transform so
// selection shapes can be mapped back into data space for the O(N) test.
type Bounds = { minX: number; maxX: number; minY: number; maxY: number };

function pointInPolygon(px: number, py: number, poly: { x: number; y: number }[]) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
        const xi = poly[i].x, yi = poly[i].y, xj = poly[j].x, yj = poly[j].y;
        if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside;
    }
    return inside;
}
function screenToData(px: number, py: number, b: Bounds, size: { w: number; h: number }, t: Transform) {
    const dataW = b.maxX - b.minX, dataH = b.maxY - b.minY;
    const s  = (dataW / dataH > size.w / size.h) ? size.w / dataW : size.h / dataH;
    const sx = s * 2 / size.w, sy = s * 2 / size.h;
    const cfx = 2 * px / size.w - 1;
    const cfy = 1 - 2 * py / size.h;
    const cx  = (cfx - 2 * t.x / size.w) / t.scale;
    const cy  = (cfy + 2 * t.y / size.h) / t.scale;
    return { x: (t.invertX ? -1 : 1) * cx / sx + b.minX + dataW / 2, y: (t.invertY ? -1 : 1) * cy / sy + b.minY + dataH / 2 };
}


const INIT_TRANSFORM: Transform = { x: 0, y: 0, scale: 1 };

export default function Viewer2d(_params: Viewer2dParams) {
    const containerRef = useRef<HTMLDivElement>(null);
    const [size, setSize]                   = useState({ w: 0, h: 0 });
    const [transform, setTransform]         = useState<Transform>(INIT_TRANSFORM);
    const [mode, setMode]                   = useState<CanvasMode>("pan");
    const [selectionMask, setSelectionMask] = useState<Float32Array | null>(null);
    const [parsedData, setParsedData]       = useState<ParsedCsvData | null>(null);

    // Refs so callbacks always see current values without being in their deps.
    const sizeRef = useRef(size); sizeRef.current = size;
    const transformRef = useRef(transform); transformRef.current = transform;
    const selectionMaskRef = useRef(selectionMask); selectionMaskRef.current = selectionMask;

    useEffect(() => {
        const ro = new ResizeObserver(([entry]) => {
            const { width, height } = entry.contentRect;
            setSize({ w: width, h: height });
        });
        ro.observe(containerRef.current!);
        return () => ro.disconnect();
    }, []);

    // Keyboard shortcuts: 1=pan, 2=rect select, 3=polygon, Escape=pan.
    useEffect(() => {
        const handler = (e: KeyboardEvent) => {
            if (e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement) return;
            if (e.key === "1") setMode("pan");
            else if (e.key === "2") setMode("select");
            else if (e.key === "3") setMode("polygon");
            else if (e.key === "Escape") setMode("pan");
        };
        window.addEventListener("keydown", handler);
        return () => window.removeEventListener("keydown", handler);
    }, []);

    const { data, dataLines, csvHeader, bounds, numClusters } = useMemo(() => {
        if (!parsedData) return { data: [] as _Cell[], dataLines: [] as string[], csvHeader: "", bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 }, numClusters: 1 };
        const { data, dataLines, csvHeader } = parsedData;
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity, n = 0;
        for (const d of data) {
            if (d.x < minX) minX = d.x; if (d.x > maxX) maxX = d.x;
            if (d.y < minY) minY = d.y; if (d.y > maxY) maxY = d.y;
            if (d.cluster > n) n = d.cluster;
        }
        return { data, dataLines, csvHeader, bounds: { minX, maxX, minY, maxY }, numClusters: n + 1 };
    }, [parsedData]);

    // byClusterEncoder — all cells colored by cluster hue. Ignores a_selected.
    const byClusterEncoder = useMemo((): ColorEncoder<_Cell> => ({
        attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
        uniforms:   [{ name: "u_numClusters", type: "float", value: numClusters }],
        colorGlsl:  `return vec4(hue2rgb(a_cluster / u_numClusters), 0.9);`,
    }), [numClusters]);

    // bySelectedEncoder — selected cells keep cluster hue, non-selected are greyed.
    // a_selected is a built-in varying provided by the selection VBO (always available).
    const bySelectedEncoder = useMemo((): ColorEncoder<_Cell> => ({
        attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
        uniforms:   [{ name: "u_numClusters", type: "float", value: numClusters }],
        colorGlsl: `
            if (a_selected < 0.5) return vec4(0.35, 0.35, 0.35, 0.35);
            return vec4(hue2rgb(a_cluster / u_numClusters), 0.9);
        `,
    }), [numClusters]);

    // O(N) membership test: convert screen rect → data coords, scan data array.
    const onSelect = useCallback((rect: { x1: number; y1: number; x2: number; y2: number }) => {
        const p1 = screenToData(rect.x1, rect.y1, bounds, sizeRef.current, transformRef.current);
        const p2 = screenToData(rect.x2, rect.y2, bounds, sizeRef.current, transformRef.current);
        const xMin = Math.min(p1.x, p2.x), xMax = Math.max(p1.x, p2.x);
        const yMin = Math.min(p1.y, p2.y), yMax = Math.max(p1.y, p2.y);
        const mask = new Float32Array(data.length);
        for (let i = 0; i < data.length; i++)
            mask[i] = data[i].x >= xMin && data[i].x <= xMax && data[i].y >= yMin && data[i].y <= yMax ? 1 : 0;
        setSelectionMask(mask);
    }, [bounds, data]);

    // O(N) membership test: convert screen polygon → data coords, ray-cast each cell.
    const onSelectPolygon = useCallback((screenPoints: { x: number; y: number }[], additive: boolean) => {
        const poly = screenPoints.map(p => screenToData(p.x, p.y, bounds, sizeRef.current, transformRef.current));
        const mask = new Float32Array(data.length);
        for (let i = 0; i < data.length; i++)
            mask[i] = pointInPolygon(data[i].x, data[i].y, poly) ? 1 : 0;
        if (additive && selectionMaskRef.current)
            for (let i = 0; i < mask.length; i++)
                if (selectionMaskRef.current[i]) mask[i] = 1;
        setSelectionMask(mask);
    }, [bounds, data]);

    const exportCsv = useCallback(() => {
        const mask     = selectionMaskRef.current;
        const baseName = parsedData?.fileName.replace(/\.csv$/i, "") ?? "cells";
        const csv      = csvHeader + ",selection_mask\n" + dataLines.map((line, i) => `${line},${mask ? mask[i] : 0}`).join("\n");
        const url      = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
        const a        = Object.assign(document.createElement("a"), { href: url, download: `${baseName}_annotated.csv` });
        a.click();
        URL.revokeObjectURL(url);
    }, [csvHeader, dataLines, parsedData]);

    const selectedCount = selectionMask ? selectionMask.reduce((n, v) => n + v, 0) : 0;

    const handleLoad = useCallback((parsed: ParsedCsvData) => {
        setParsedData(parsed);
        setSelectionMask(parsed.initialMask);
        setTransform(INIT_TRANSFORM);
    }, []);

    return (
        <div ref={containerRef} style={{ position: "relative", width: "100%", height: "100vh" }}>
            {!parsedData && <CsvUploadDialog onLoad={handleLoad} />}
            <UnderlyingCanvas
                data={data}
                xAccessor={d => d.x}
                yAccessor={d => d.y}
                colorEncoder={selectionMask ? bySelectedEncoder : byClusterEncoder}
                transform={transform}
                size={size}
                selectionMask={selectionMask}
            />
            <OverlyingCanvas
                size={size}
                mode={mode}
                transform={transform}
                onTransform={setTransform}
                onSelect={onSelect}
                onSelectPolygon={onSelectPolygon}
            />
            <div style={{ position: "absolute", top: 12, left: 12, display: "flex", gap: 6, alignItems: "center", color: "white", fontFamily: "sans-serif", fontSize: 14 }}>
                {(["pan", "select", "polygon"] as CanvasMode[]).map((m, i) => (
                    <button key={m} onClick={() => setMode(m)}
                        style={{ opacity: mode === m ? 1 : 0.5, fontWeight: mode === m ? "bold" : "normal" }}>
                        {m === "pan" ? "Pan" : m === "select" ? "Rect" : "Polygon"} [{i + 1}]
                    </button>
                ))}
                <span style={{ color: "rgba(255,255,255,0.4)", margin: "0 4px" }}>|</span>
                <button onClick={() => setTransform(t => ({ ...t, invertX: !t.invertX }))}
                    style={{ opacity: transform.invertX ? 1 : 0.5, fontWeight: transform.invertX ? "bold" : "normal" }}>
                    Inv X
                </button>
                <button onClick={() => setTransform(t => ({ ...t, invertY: !t.invertY }))}
                    style={{ opacity: transform.invertY ? 1 : 0.5, fontWeight: transform.invertY ? "bold" : "normal" }}>
                    Inv Y
                </button>
                <span style={{ color: "rgba(255,255,255,0.4)", margin: "0 4px" }}>|</span>
                {selectionMask && (
                    <span style={{ color: "white", fontSize: 12 }}>{selectedCount | 0} selected</span>
                )}
                {selectionMask && (
                    <button onClick={() => setSelectionMask(null)}>Clear</button>
                )}
                <button onClick={exportCsv}>Export CSV</button>
            </div>
        </div>
    );
}