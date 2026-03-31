"use client";

import { ScatterCanvas } from "../brush-2d-v2/scatter-canvas";
import HorizontalSplit from "../horizontal-split";

export default function Page() {
    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            <HorizontalSplit>
                <div style={{ background: "lightblue" }}>
                    <ScatterCanvas
                        key={lcSubset}
                        data={pointsDataSubset}
                        xAccessor={}
                        yAccessor={}
                        colorEncoder={byClusterSelectionEncoderV2}
                        size={lcSize}
                        transform={ lcTransform }
                        setTransform={setLcTransform}
                        onReady={() => {}}
                        polygonMask={polygonMask}
                    ></ScatterCanvas>
                </div>
                <div style={{ background: "lightcoral" }}>Hello world right</div>
            </HorizontalSplit>
        </div>
    );
}