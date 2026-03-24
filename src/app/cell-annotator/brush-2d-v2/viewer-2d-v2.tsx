"use client";

import { useEffect, useRef, useState } from "react";
import { ColorEncoder, Transform } from "../underlying-canvas";
import PolygonManagerCanvas from "./polygon-manager-canvas";
import type { PolygonManagerHandle } from "./polygon-manager-canvas";
import { OverlyingCanvasV2 } from "./overlying-canvas-v2";
import { ScatterCanvas } from "./scatter-canvas";

interface CellData {
    id: string;
    x: number;
    y: number;
    cluster: number;
}

const ALL_PATIENTS = [
    "LARC_A_S19_31776B1", "LARC_A_S20_12903C1", "LARC_A_S21_12210A3",
    "LARC_A_S21_1580B1", "LARC_A_S21_19043A1", "LARC_A_S22_1068A1",
    "LARC_A_S22_12362A1", "LARC_A_SP18_2074B1", "LARC_A_SP21_204B1",
    "LARC_B_S19_12126B1", "LARC_B_S19_25142A1", "LARC_B_S20_19504A1",
    "LARC_B_S20_23887A1", "LARC_B_S20_25858B1", "LARC_B_S21_12126A1",
    "LARC_B_S21_8615C1", "LARC_B_SP19_4321C1", "LARC_B_SP22_1201A1",
    "LARC_B_UNIDENTIFIED",
];

async function loadPatientCsv(patientId: string): Promise<CellData[]> {
    const res = await fetch(`/cell-annotator-data/${patientId}_annotated.csv`);
    const text = await res.text();
    return text.trim().split("\n").slice(1).map(line => {
        const [id, x, y, cluster] = line.split(",");
        return { id, x: +x, y: +y, cluster: +cluster };
    });
}

const byClusterEncoder: ColorEncoder<CellData> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
};

const byClusterSelectionEncoder: ColorEncoder<CellData> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `
        if (a_polygon == 0.0) return vec4(0.35, 0.35, 0.35, 1.0);
        return vec4(hue2rgb(a_cluster / 23.0), 1.0);
    `,
};

const INIT_TRANSFORM: Transform = { x: 0, y: 0, scale: 1 };

function fitTransform(data: CellData[], w: number, h: number): Transform {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const d of data) {
        if (d.x < minX) minX = d.x; if (d.x > maxX) maxX = d.x;
        if (d.y < minY) minY = d.y; if (d.y > maxY) maxY = d.y;
    }
    const dataW = maxX - minX, dataH = maxY - minY;
    const scale = dataW / dataH > w / h ? w / dataW : h / dataH;
    const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
    return { x: w / 2 - cx * scale, y: h / 2 + cy * scale, scale };
}

function SubsetSelector({ patients, selected, onSelect }: {
    patients: string[];
    selected: string;
    onSelect: (p: string) => void;
}) {
    return (
        <div className="thin-scrollbar" style={{
            position: "absolute", bottom: 0, left: "50%", transform: "translateX(-50%)",
            display: "flex", gap: 6, padding: "6px 10px",
            background: "rgba(0,0,0,0.55)", backdropFilter: "blur(6px)",
            borderRadius: 10, overflowX: "auto", maxWidth: "100%",
            boxShadow: "0 2px 12px rgba(0,0,0,0.4)",
        }}>
            {patients.map(p => (
                <button key={p} onClick={() => onSelect(p)} style={{
                    flexShrink: 0, padding: "4px 12px", borderRadius: 6, border: "none",
                    cursor: "pointer", fontSize: 12, fontFamily: "monospace",
                    background: p === selected ? "rgba(255,255,255,0.9)" : "rgba(255,255,255,0.12)",
                    color: p === selected ? "#111" : "#ddd",
                    fontWeight: p === selected ? 700 : 400,
                }}>
                    {p}
                </button>
            ))}
        </div>
    );
}

export default function Viewer2d() {
    const containerRef = useRef<HTMLDivElement>(null);
    const polygonManagerRef = useRef<PolygonManagerHandle | null>(null);
    const [size, setSize]           = useState({ w: 0, h: 0 });
    const [transform, setTransform] = useState<Transform>(INIT_TRANSFORM);
    const [subset, setSubset]             = useState(ALL_PATIENTS[0]);
    const [dataSubset, setDataSubset]     = useState<CellData[]>([]);
    const [loadedSubset, setLoadedSubset] = useState("");
    const [readySubset, setReadySubset]   = useState("");
    const [rightWidth, setRightWidth] = useState(500);
    const dragRef = useRef<{ startX: number; startW: number } | null>(null);

    const sizeRef      = useRef(size);      sizeRef.current      = size;
    const transformRef = useRef(transform); transformRef.current = transform;
    const fittedDataRef = useRef<CellData[] | null>(null);

    useEffect(() => {
        loadPatientCsv(subset).then(data => {
            setDataSubset(data);
            setLoadedSubset(subset);
        });
    }, [subset]);

    useEffect(() => {
        if (loadedSubset !== subset || dataSubset.length === 0 || size.w === 0 || fittedDataRef.current === dataSubset) return;
        fittedDataRef.current = dataSubset;
        setTransform(fitTransform(dataSubset, size.w, size.h));
    }, [dataSubset, loadedSubset, size, subset]);

    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    const onSelect = (patientId: string) => {
        setSubset(patientId);
    }

    const onDividerMouseDown = (e: React.MouseEvent) => {
        e.preventDefault();
        dragRef.current = { startX: e.clientX, startW: rightWidth };
        const onMove = (ev: MouseEvent) => {
            const delta = dragRef.current!.startX - ev.clientX;
            setRightWidth(Math.max(200, Math.min(900, dragRef.current!.startW + delta)));
        };
        const onUp = () => {
            dragRef.current = null;
            window.removeEventListener("mousemove", onMove);
            window.removeEventListener("mouseup", onUp);
        };
        window.addEventListener("mousemove", onMove);
        window.addEventListener("mouseup", onUp);
    };

    const scatterData = loadedSubset === subset ? dataSubset : [];
    const isLoading = readySubset !== subset;

    return (
        <div style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
            <div ref={containerRef} style={{ position: "relative", flex: "1 1 auto", background: "#111" }}>
                <ScatterCanvas
                    key={subset}  // Force remount when subset changes, to reset WebGL state
                    data={scatterData}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={byClusterEncoder}
                    size={size}
                    transform={transform}
                    onReady={() => setReadySubset(subset)}
                />
                <PolygonManagerCanvas
                    handleRef={polygonManagerRef}
                    size={size}
                    transform={transform}
                    onTransform={setTransform}
                    subset={subset}
                />
                <OverlyingCanvasV2
                    size={size}
                    mode="brush"
                    transform={transform}
                    onTransform={setTransform}
                    onBrush={(x, y) => polygonManagerRef.current?.onBrushClick(x, y)}
                    onBrushMove={(x, y) => polygonManagerRef.current?.onBrushMove(x, y)}
                    onCursorMove={(x, y) => polygonManagerRef.current?.setCursorPos(x, y)}
                    onCursorLeave={() => polygonManagerRef.current?.clearCursor()}
                    onBrushResize={delta => polygonManagerRef.current?.adjustBrushRadius(delta)}
                />
                <div style={{
                    position: "absolute", inset: 0, zIndex: 50,
                    background: "rgba(17,17,17,0.5)",
                    backdropFilter: "blur(14px)",
                    WebkitBackdropFilter: "blur(14px)",
                    border: "1px solid rgba(255,255,255,0.06)",
                    transition: isLoading ? "none" : "opacity 0.35s ease",
                    opacity: isLoading ? 1 : 0,
                    pointerEvents: isLoading ? "all" : "none",
                }} />
            </div>
            <div
                onMouseDown={onDividerMouseDown}
                style={{
                    width: 5, flexShrink: 0, cursor: "col-resize",
                    background: "rgba(255,255,255,0.08)",
                }}
                onMouseEnter={e => (e.currentTarget.style.background = "rgba(255,255,255,0.22)")}
                onMouseLeave={e => (e.currentTarget.style.background = "rgba(255,255,255,0.08)")}
            />
            <div style={{ flex: `0 0 ${rightWidth}px`, background: "#222", color: "#fff", padding: 16, position: 'relative'}}>
                {/*<p style={{ fontFamily: "monospace", fontSize: 12 }}>{subset}</p>
                <p style={{ fontFamily: "monospace", fontSize: 12 }}>{isLoading ? "Loading…" : `Points: ${dataSubset.length}`}</p>
                */}
                
                <ScatterCanvas
                    data={scatterData}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={byClusterSelectionEncoder}
                    size={{ w: rightWidth, h: size.h }}
                    transform={transform}
                ></ScatterCanvas>
            </div>
            <SubsetSelector patients={ALL_PATIENTS} selected={subset} onSelect={onSelect} />
        </div>
    );
}
