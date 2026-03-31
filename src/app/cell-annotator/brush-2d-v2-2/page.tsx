"use client";

import { useEffect, useRef, useState } from "react";
import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import CsvUploadDialogV2, { CellDataV2 } from "../brush-2d-v2/csv-upload-dialog-v2";
import { ColorEncoder, Transform } from "../underlying-canvas";
import { byClusterSelectionEncoder } from "../brush-2d-v2/viewer-2d-v2.1";

const colorEncoder: ColorEncoder<CellDataV2> = {
    attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
    uniforms: [],
    colorGlsl: `return vec4(hue2rgb(a_cluster / 23.0), 1.0);`,
};

export default function Page() {
    const containerRef = useRef<HTMLDivElement>(null);
    const [size, setSize] = useState({ w: 0, h: 0 });
    const [lcTransform, setLcTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [pointsData, setPointsData] = useState<CellDataV2[]>([]);

    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            {pointsData.length === 0 && <CsvUploadDialogV2 onLoad={setPointsData} />}
            <div ref={containerRef} style={{ width: "100%", height: "100%" }}>
                <ScatterCanvas
                    data={pointsData}
                    xAccessor={d => d.x}
                    yAccessor={d => d.y}
                    colorEncoder={byClusterSelectionEncoder}
                    size={size}
                    transform={lcTransform}
                    setTransform={setLcTransform}
                    onReady={() => {}}
                />
            </div>
        </div>
    );
}
