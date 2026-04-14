"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import CsvUploadDialogV2_1, { CellDataV2_1 } from "../brush-2d-v2/csv-upload-dialog-v2-1";
import { ColorEncoder, Transform } from "../underlying-canvas";
import { mat4mul, perspectiveMat, viewMat } from "../brush-3d/scatter-3d";
import SubsetSelectorV2_1 from "../brush-2d-v2/subset-selector-v2.1";
import ClusterSelectorV2_1 from "../brush-2d-v2/cluster-selector-v2-1";
import { OverlyingCanvasV2 } from "../brush-2d-v2/overlying-canvas-v2";
import { OverlyingCanvasV3 } from "../brush-2d-v2/overlying-canvas-v3";
import DotSizeSelector from "./dot-size-selector";
import HorizontalSplit from "../horizontal-split";
import { Scatter3DV2 } from "./scatter-3d-v2";

const colorEncoder: ColorEncoder<CellDataV2_1> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.clusterIdx }],
    uniforms: [],
    colorGlsl: `
        if (a_polygon == 0.0) return vec4(0.35, 0.35, 0.35, 1.0);
        return vec4(hue2rgb(a_cluster / 15.0), 1.0);
    `,
};

const byClusterColorEncoder: ColorEncoder<CellDataV2_1> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.clusterIdx }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 15.0), 1.0);`,
};

const byClusterAndSelectionColorEncoder: ColorEncoder<CellDataV2_1> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.clusterIdx }],
    uniforms: [],
    colorGlsl: `
        if (a_polygon == 0.0) return vec4(0.35, 0.35, 0.35, 1.0);
        return vec4(hue2rgb(a_cluster / 15.0), 1.0);
    `,
};

type _CellData = CellDataV2_1;

export default function Page() {
    const lcRef = useRef<HTMLDivElement>(null);
    const [lcSize, setLcSize] = useState({ w: 0, h: 0 });
    const [lcTransform, setLcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [pointsData, setPointsData] = useState<CellDataV2_1[]>([]);
    const [groupList, setGroupList] = useState<string[]>([]); // patient groups for subset selector

    const [lcSubset, setLcSubset] = useState<string | null>(null);
    const [pointsDataSubset, setPointsDataSubset] = useState<_CellData[]>([]);
    const [selectionMask, setSelectionMask] = useState<Float32Array | null>(null);
    const [lcDotSize, setLcDotSize] = useState(2);


    const rcRef = useRef<HTMLDivElement>(null);
    const [rcSize, setRcSize] = useState({ w: 0, h: 0 });
    const [rcMatrix, setRcMatrix] = useState(() => new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]));

    const rcAspect = rcSize.w > 0 ? rcSize.w / rcSize.h : 1;
    const rcFullMatrix = mat4mul(perspectiveMat(60 * Math.PI / 180, rcAspect, 0.1, 100.0), mat4mul(viewMat(3.0), rcMatrix));

    // Mirror the normalization done inside Scatter3DV2 so OverlyingCanvasV3 works in the same space.
    const rcNormalizedPoints = useMemo(() => {
        let minX = Infinity, maxX = -Infinity;
        let minY = Infinity, maxY = -Infinity;
        let minZ = Infinity, maxZ = -Infinity;
        for (const d of pointsDataSubset) {
            if (d.umap_1 < minX) minX = d.umap_1; if (d.umap_1 > maxX) maxX = d.umap_1;
            if (d.umap_2 < minY) minY = d.umap_2; if (d.umap_2 > maxY) maxY = d.umap_2;
            const z = d.umap_3 ?? 0;
            if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
        }
        const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ) / 1.6 || 1;
        const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2, cz = (minZ + maxZ) / 2;
        return pointsDataSubset.map(d => ({ x: (d.umap_1 - cx) / span, y: (d.umap_2 - cy) / span, z: ((d.umap_3 ?? 0) - cz) / span }));
    }, [pointsDataSubset]);

    const handleLoad = (data: CellDataV2_1[]) => {
        console.log("Data", data);
        setPointsData(data);
        const patients = Array.from(new Set(data.map(d => d.SampleId)));
        setGroupList(patients);
        setLcSubset(patients[0]);
    };

    useEffect(() => {
        if (lcSubset) setPointsDataSubset(pointsData.filter(d => d.SampleId === lcSubset));
    }, [lcSubset, pointsData]);

    useEffect(() => {
        const el = lcRef.current!;
        const ro = new ResizeObserver(() => setLcSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);

        const rEl = rcRef.current!;
        const rRo = new ResizeObserver(() => setRcSize({ w: rEl.clientWidth, h: rEl.clientHeight }));
        rRo.observe(rEl);
        return () => {
            ro.disconnect();
            rRo.disconnect();
        };
    }, []);

    // Debugging only
    useEffect(() => {
        if (!selectionMask) return;

        const selectedPoints = pointsDataSubset.filter((_, i) => selectionMask[i] === 1);
        const selectedIds = new Set(selectedPoints.map(p => p.id));

        const text = Array.from(selectedIds).join("\n");

        navigator.clipboard.writeText(text)
            .then(() => console.log("Copied to clipboard"))
            .catch(err => console.error("Copy failed:", err));

    }, [pointsDataSubset, selectionMask]);

    // FOR DEBUGGIN PURPOSE ONLY
    /*useEffect(() => {
        fetch("/cell-annotator-data-v2/umap_clusters_3d.csv")
            .then(r => r.text())
            .then(text => {
                const lines = text.split(/\r?\n/);
                const headers = lines[0].split(",").map(h => h.replace(/"/g, ""));
                const [xIdx, yIdx, u1Idx, u2Idx, u3Idx, cIdx, sIdx, idIdx] = ["x","y","umap_1","umap_2","umap_3","cluster","SampleId","cell"].map(n => headers.indexOf(n));
                const clusterMap = new Map<string, number>();
                const clusterOrder: string[] = [];
                const data = lines.slice(1).filter(Boolean).map((line, i) => {
                    const cols = line.split(",").map(c => c.replace(/"/g, ""));
                    const cluster = cols[cIdx];
                    if (!clusterMap.has(cluster)) { clusterMap.set(cluster, clusterOrder.length); clusterOrder.push(cluster); }
                    return { id: idIdx >= 0 ? cols[idIdx] : String(i), x: +cols[xIdx], y: +cols[yIdx], umap_1: +cols[u1Idx], umap_2: +cols[u2Idx], umap_3: +cols[u3Idx], cluster, clusterIdx: clusterMap.get(cluster)!, SampleId: cols[sIdx] };
                });
                handleLoad(data);
            });
    }, []);*/


    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            {pointsData.length === 0 && <CsvUploadDialogV2_1 onLoad={handleLoad} />}
            <HorizontalSplit>
                <div ref={lcRef} style={{ height: "100%" }}>
                    <ClusterSelectorV2_1
                        data={pointsDataSubset}
                        colorEncoder={colorEncoder}
                        onMaskChange={setSelectionMask}
                    />
                    <ScatterCanvas
                        key={lcSubset}
                        data={pointsDataSubset}
                        xAccessor={d => d.x}
                        yAccessor={d => d.y}
                        colorEncoder={byClusterAndSelectionColorEncoder}
                        size={lcSize}
                        transform={lcTransform}
                        setTransform={setLcTransform}
                        onReady={() => {}}
                        polygonMask={selectionMask}
                        dotSize={lcDotSize}
                    />
                    <OverlyingCanvasV2
                        size={lcSize}
                        mode="pan"
                        transform={lcTransform}
                        onTransform={setLcTransform}
                        onBrush={() => {}}
                    ></OverlyingCanvasV2>
                    <SubsetSelectorV2_1
                        subsets={groupList}
                        subset={lcSubset}
                        onSubsetChange={setLcSubset}
                    ></SubsetSelectorV2_1>
                    <DotSizeSelector dotSize={lcDotSize} onDotSizeChange={setLcDotSize} />
                </div>
                <div ref={rcRef} style={{ height: "100%" }}>
                    <Scatter3DV2
                        data={pointsDataSubset}
                        xAccessor={d => d.umap_1}
                        yAccessor={d => d.umap_2}
                        zAccessor={d => d.umap_3 ?? 0}
                        colorEncoder={byClusterAndSelectionColorEncoder}
                        polygonMask={selectionMask}
                        size={rcSize}
                        matrix={rcFullMatrix}
                        onReady={() => {}}
                    ></Scatter3DV2>
                    <OverlyingCanvasV3
                        size={rcSize}
                        matrix={rcFullMatrix}
                        modelMatrix={rcMatrix}
                        points={rcNormalizedPoints}
                        onBrush={indices => {
                            const mask = new Float32Array(pointsDataSubset.length);
                            indices.forEach(i => { mask[i] = 1.0; });
                            setSelectionMask(mask);
                        }}
                        onModelMatrixChange={m => setRcMatrix(m as Float32Array<ArrayBuffer>)}
                    />
                    <div style={{ position: "absolute", bottom: 16, right: 16, display: "flex", gap: 4, zIndex: 10, color: 'white' }}>
                        <button onClick={() => setRcMatrix(m => mat4mul(new Float32Array([1.2,0,0,0, 0,1.2,0,0, 0,0,1.2,0, 0,0,0,1]), m) as Float32Array<ArrayBuffer>)}>[+]</button>
                        <button onClick={() => setRcMatrix(m => mat4mul(new Float32Array([1/1.2,0,0,0, 0,1/1.2,0,0, 0,0,1/1.2,0, 0,0,0,1]), m) as Float32Array<ArrayBuffer>)}>[-]</button>
                    </div>
                </div>
            </HorizontalSplit>
        </div>
    );
}
