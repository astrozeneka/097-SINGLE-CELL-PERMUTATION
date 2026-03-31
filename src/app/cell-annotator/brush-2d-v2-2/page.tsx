"use client";

import { useEffect, useRef, useState } from "react";
import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import CsvUploadDialogV2_1, { CellDataV2_1 } from "../brush-2d-v2/csv-upload-dialog-v2-1";
import { ColorEncoder, Transform } from "../underlying-canvas";
import SubsetSelectorV2_1 from "../brush-2d-v2/subset-selector-v2.1";

const colorEncoder: ColorEncoder<CellDataV2_1> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.clusterIdx }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
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

    const handleLoad = (data: CellDataV2_1[]) => {
        console.log("===data", data);
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
                />
                <SubsetSelectorV2_1
                    subsets={groupList}
                    subset={lcSubset}
                    onSubsetChange={setLcSubset}
                ></SubsetSelectorV2_1>
            </div>
        </div>
    );
}
