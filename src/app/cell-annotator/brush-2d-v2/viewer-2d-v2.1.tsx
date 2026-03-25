"use client";

import { useEffect, useRef, useState } from "react";
import { OverlyingCanvasV2 } from "./overlying-canvas-v2";
import PolygonManagerCanvas, { PolygonManagerHandle } from "./polygon-manager-canvas";
import { ScatterCanvas } from "./scatter-canvas";
import { ColorEncoder, Transform } from "../underlying-canvas";

interface _CellData {
    id: string;
    x: number;
    y: number;
    umap_1: number;
    umap_2: number;
    cluster: number;
    SampleId: string;
}

async function loadAllPatientsCsv(): Promise<_CellData[]> {
    const res = await fetch(`/cell-annotator-data-v2/cell_coordinates_with_sample.csv`, { cache: "no-store" });
    const text = await res.text();
    return text.trim().split("\n").slice(1).map(line => {
        const [id, x, y, umap_1, umap_2, cluster, SampleId] = line.split(",");
        return { id, x: +x, y: +y, umap_1: +umap_1, umap_2: +umap_2, cluster: +cluster, SampleId };

    });
}

const byClusterEncoder: ColorEncoder<_CellData> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
};

export default function Viewer2dPCA_viewer() {
    const lcRef = useRef<HTMLDivElement>(null);
    const rcRef = useRef<HTMLDivElement>(null);

    const polygonManagerRef = useRef<PolygonManagerHandle | null>(null);
    const [lcSize, setLcSize]           = useState({ w: 0, h: 0 });
    const [rcSize, setRcSize]           = useState({ w: 0, h: 0 });
    const [rcTransform, setRcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [subset, setSubset] = useState("default");

    const [rightWidth, setRightWidth] = useState(500);
    const [polygons, setPolygons] = useState<{ verts: { x: number; y: number }[] }[]>([]);
    const dragRef = useRef<{ startX: number; startW: number } | null>(null);

    // Newly introduced int he v2.1
    const [pointsData, setPointsData] = useState<_CellData[]>([]);

    useEffect(() => {
        loadAllPatientsCsv().then(data => {
            setPointsData(data);
            // setLoadedSubset("all"); 
        });
    }, []);

    useEffect(() => {
        // Fit, later, this should be handled by the scatter canvas itself
        // Not here
    })

    useEffect(() => {
        const lc = lcRef.current!;
        const lcRo = new ResizeObserver(() => setLcSize({ w: lc.clientWidth, h: lc.clientHeight }));
        lcRo.observe(lc);


        const rc = rcRef.current!;
        const rcRo = new ResizeObserver(() => setRcSize({ w: rc.clientWidth, h: rc.clientHeight }));
        rcRo.observe(rc);

        return () => rcRo.disconnect();
    }, []);


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

    useEffect(()=>{
        console.log("T", rcTransform);
    }, [rcTransform]);

    return (
        
        <div style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
            <div ref={lcRef} style={{ position: "relative", flex: "1 1 auto", background: "#111" }}>


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
            <div
                ref={rcRef}
                style={{ flex: `0 0 ${rightWidth}px`, background: "#222", color: "#fff", padding: 16, position: 'relative'}}
            >
                <ScatterCanvas
                    key={subset}
                    data={pointsData}
                    xAccessor={d => d.umap_1}
                    yAccessor={d => d.umap_2}
                    colorEncoder={byClusterEncoder}
                    size={rcSize}
                    transform={rcTransform}
                    onReady={() => {}}
                ></ScatterCanvas>
                <PolygonManagerCanvas
                    handleRef={polygonManagerRef}
                    size={rcSize}
                    transform={rcTransform}
                    onTransform={setRcTransform}
                    subset={subset}
                    onPolygonsChange={setPolygons}
                ></PolygonManagerCanvas>
                <OverlyingCanvasV2
                    size={rcSize}
                    mode="brush"
                    transform={rcTransform}
                    onTransform={setRcTransform}
                    onBrush={() => {}}
                    
                ></OverlyingCanvasV2>
            </div>
        </div>
    );
}