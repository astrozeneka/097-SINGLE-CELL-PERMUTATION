"use client";

import { useEffect, useRef, useState, useMemo } from "react";
import { OverlyingCanvasV2 } from "./overlying-canvas-v2";
import PolygonManagerCanvas, { PolygonManagerHandle, pointInPolygon } from "./polygon-manager-canvas";
import { ScatterCanvas } from "./scatter-canvas";
import { ColorEncoder, Transform } from "../underlying-canvas";
import SubsetSelectorV2_1 from "./subset-selector-v2.1";
import CsvUploadDialogV2, { CellDataV2 } from "./csv-upload-dialog-v2";
import ClusterSelector from "./cluster-selector";

type _CellData = CellDataV2;

const byClusterEncoder: ColorEncoder<_CellData> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
};

export const byClusterSelectionEncoder: ColorEncoder<_CellData> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `
        if (a_polygon >= 0.0) return vec4(0.35, 0.35, 0.35, 1.0);
        return vec4(hue2rgb(a_cluster / 23.0), 1.0);
    `,
};

export default function Viewer2dPCA_viewer() {
    const lcRef = useRef<HTMLDivElement>(null);
    const rcRef = useRef<HTMLDivElement>(null);

    const polygonManagerRef = useRef<PolygonManagerHandle | null>(null);
    const [lcSize, setLcSize]           = useState({ w: 0, h: 0 });
    const [rcSize, setRcSize]           = useState({ w: 0, h: 0 });
    const [lcTransform, setLcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [rcTransform, setRcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });

    const [rightWidth, setRightWidth] = useState(500);
    const [polygons, setPolygons] = useState<{ verts: { x: number; y: number }[] }[]>([]);
    const dragRef = useRef<{ startX: number; startW: number } | null>(null);

    // Newly introduced int he v2.1
    const [pointsData, setPointsData] = useState<_CellData[]>([]);

    const [patients, setPatients] = useState<string[]>([]);

    // Left container
    const [lcSubset, setLcSubset] = useState<string | null>(null);
    const [pointsDataSubset, setPointsDataSubset] = useState<_CellData[]>([]);


    const handleLoad = (data: _CellData[]) => {
        setPointsData(data);
        const patients = Array.from(new Set(data.map(d => d.SampleId)));
        setPatients(patients);
        setLcSubset(patients[0]);
    };

    useEffect(() => {
        if (lcSubset) {
            setPointsDataSubset(pointsData.filter(d => d.SampleId === lcSubset));
        }
    }, [lcSubset, pointsData]);

    useEffect(() => {
        const lc = lcRef.current!;
        const lcRo = new ResizeObserver(() => setLcSize({ w: lc.clientWidth, h: lc.clientHeight }));
        lcRo.observe(lc);


        const rc = rcRef.current!;
        const rcRo = new ResizeObserver(() => setRcSize({ w: rc.clientWidth, h: rc.clientHeight }));
        rcRo.observe(rc);

        return () => rcRo.disconnect();
    }, []);

    const polygonMask = useMemo(() => {
        if (polygons.length === 0) return null;
        // AN idea lcase is to put this inside the scatter canvas ?? (maybe not)
        const mask = new Float32Array(pointsDataSubset.length);
        for (let i = 0; i < pointsDataSubset.length; i++) {
            const { umap_1, umap_2 } = pointsDataSubset[i];
            if (polygons.some(p => pointInPolygon(p.verts, umap_1, umap_2))) mask[i] = 1.0;
        }
        console.log("Mask", mask)
        return mask;
    }, [pointsDataSubset, polygons]);


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

    return (
        
        <div style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
            {pointsData.length === 0 && <CsvUploadDialogV2 onLoad={handleLoad} />}
            <div ref={lcRef} style={{ position: "relative", flex: "1 1 auto", background: "#111" }}>
                
                
                <ScatterCanvas
                    key={lcSubset}
                    data={pointsDataSubset}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={byClusterSelectionEncoder}
                    size={lcSize}
                    transform={ lcTransform }
                    setTransform={setLcTransform}
                    onReady={() => {}}
                    polygonMask={polygonMask}
                ></ScatterCanvas>
                <OverlyingCanvasV2
                    size={lcSize}
                    mode="pan"
                    transform={lcTransform}
                    onTransform={setLcTransform}
                    onBrush={() => {}}
                ></OverlyingCanvasV2>
                <SubsetSelectorV2_1
                    subsets={patients}
                    subset={lcSubset}
                    onSubsetChange={setLcSubset}
                ></SubsetSelectorV2_1>
                
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
                <ClusterSelector
                    data={pointsData}
                    colorEncoder={byClusterEncoder}
                ></ClusterSelector>
                <ScatterCanvas
                    key="default"
                    data={pointsData}
                    xAccessor={d => d.umap_1}
                    yAccessor={d => d.umap_2}
                    colorEncoder={byClusterEncoder}
                    size={rcSize}
                    transform={rcTransform}
                    setTransform={setRcTransform}
                    onReady={() => {}}
                ></ScatterCanvas>
                <PolygonManagerCanvas
                    handleRef={polygonManagerRef}
                    size={rcSize}
                    transform={rcTransform}
                    onTransform={setRcTransform}
                    subset="default"
                    onPolygonsChange={setPolygons}
                ></PolygonManagerCanvas>
                <OverlyingCanvasV2
                    size={rcSize}
                    mode="brush"
                    transform={rcTransform}
                    onTransform={setRcTransform}
                    onBrush={(x, y) => polygonManagerRef.current?.onBrushClick(x, y)}
                    onBrushMove={(x, y) => polygonManagerRef.current?.onBrushMove(x, y)}
                    onCursorMove={(x, y) => polygonManagerRef.current?.setCursorPos(x, y)}
                    onCursorLeave={() => polygonManagerRef.current?.clearCursor()}
                    onBrushResize={delta => polygonManagerRef.current?.adjustBrushRadius(delta)}
                ></OverlyingCanvasV2>
                <div style={{ position: 'absolute', bottom: 16, left: 16, zIndex: 10, fontFamily: 'monospace', fontSize: 11 }}>
                    <button onClick={() => polygonManagerRef.current?.clearPolygons()}>Clear selection</button>
                </div>
            </div>
        </div>
    );
}