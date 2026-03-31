"use client";

import { useEffect, useRef, useState } from "react";
import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import CsvUploadDialogV2, { CellDataV2 } from "../brush-2d-v2/csv-upload-dialog-v2";
import { ColorEncoder, Transform } from "../underlying-canvas";
import { byClusterSelectionEncoder } from "../brush-2d-v2/viewer-2d-v2.1";
import SubsetSelectorV2_1 from "../brush-2d-v2/subset-selector-v2.1";

const colorEncoder: ColorEncoder<CellDataV2> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
};

type _CellData = CellDataV2;

export default function Page() {
    const containerRef = useRef<HTMLDivElement>(null);
    const [size, setSize] = useState({ w: 0, h: 0 });
    const [lcTransform, setLcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [pointsData, setPointsData] = useState<CellDataV2[]>([]);
    const [groupList, setGroupList] = useState<string[]>([]); // patient groups for subset selector

    const [lcSubset, setLcSubset] = useState<string | null>(null);
    const [pointsDataSubset, setPointsDataSubset] = useState<_CellData[]>([]);

    const handleLoad = (data: CellDataV2[]) => {
        setPointsData(data);
        const patients = Array.from(new Set(data.map(d => d.SampleId)));
        setGroupList(patients);
    };

    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            {pointsData.length === 0 && <CsvUploadDialogV2 onLoad={handleLoad} />}
            <div ref={containerRef} style={{ width: "100%", height: "100%" }}>
                <ScatterCanvas
                    key={lcSubset}
                    data={pointsDataSubset}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={byClusterSelectionEncoder}
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
