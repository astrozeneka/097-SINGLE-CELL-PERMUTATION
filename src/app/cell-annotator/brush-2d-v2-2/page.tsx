"use client";

import { useEffect, useRef, useState } from "react";
import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import CsvUploadDialogV2_1, { CellDataV2_1 } from "../brush-2d-v2/csv-upload-dialog-v2-1";
import { ColorEncoder, Transform } from "../underlying-canvas";
import SubsetSelectorV2_1 from "../brush-2d-v2/subset-selector-v2.1";
import ClusterSelectorV2_1 from "../brush-2d-v2/cluster-selector-v2-1";
import { OverlyingCanvasV2 } from "../brush-2d-v2/overlying-canvas-v2";

const colorEncoder: ColorEncoder<CellDataV2_1> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.clusterIdx }],
    uniforms: [],
    colorGlsl: `
        if (a_polygon == 0.0) return vec4(0.35, 0.35, 0.35, 1.0);
        return vec4(hue2rgb(a_cluster / 15.0), 1.0);
    `,
};

type _CellData = CellDataV2_1;

export default function Page() {
    const containerRef = useRef<HTMLDivElement>(null);
    const [size, setSize] = useState({ w: 0, h: 0 });
    const [lcTransform, setLcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [pointsData, setPointsData] = useState<CellDataV2_1[]>([]);
    const [groupList, setGroupList] = useState<string[]>([]); // patient groups for subset selector

    const [lcSubset, setLcSubset] = useState<string | null>(null);
    const [pointsDataSubset, setPointsDataSubset] = useState<_CellData[]>([]);
    const [clusterMask, setClusterMask] = useState<Float32Array | null>(null);

    const handleLoad = (data: CellDataV2_1[]) => {
        setPointsData(data);
        const patients = Array.from(new Set(data.map(d => d.SampleId)));
        setGroupList(patients);
        setLcSubset(patients[0]);
    };

    useEffect(() => {
        if (lcSubset) setPointsDataSubset(pointsData.filter(d => d.SampleId === lcSubset));
    }, [lcSubset, pointsData]);

    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            {pointsData.length === 0 && <CsvUploadDialogV2_1 onLoad={handleLoad} />}
            <div ref={containerRef} style={{ width: "100%", height: "100%" }}>
                <ClusterSelectorV2_1
                    data={pointsDataSubset}
                    colorEncoder={colorEncoder}
                    onMaskChange={setClusterMask}
                />
                <ScatterCanvas
                    key={lcSubset}
                    data={pointsDataSubset}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={colorEncoder}
                    size={size}
                    transform={lcTransform}
                    setTransform={setLcTransform}
                    onReady={() => {}}
                    polygonMask={clusterMask}
                />
                <OverlyingCanvasV2
                    size={size}
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
            </div>
        </div>
    );
}
