"use client";

import { useEffect, useRef, useState } from "react";
import { Transform } from "../utils";
import PolygonManagerCanvas from "./polygon-manager-canvas";
import type { PolygonManagerHandle } from "./polygon-manager-canvas";
import { OverlyingCanvasV2 } from "./overlying-canvas-v2";


interface Viewer2dParams {
    // TODO later
}

const INIT_TRANSFORM = { x: 0, y: 0, scale: 1 };

export default function Viewer2d(_params: Viewer2dParams) {
    const containerRef = useRef<HTMLDivElement>(null);
    const polygonManagerRef = useRef<PolygonManagerHandle | null>(null);
    const [size, setSize]                   = useState({ w: 0, h: 0 });
    const [transform, setTransform]         = useState<Transform>(INIT_TRANSFORM);
   
    // Refs so callbacks always see current values without being in their deps.
    const sizeRef = useRef(size); sizeRef.current = size;
    const transformRef = useRef(transform); transformRef.current = transform;
    // const selectionMaskRef = useRef(selectionMask); selectionMaskRef.current = selectionMask; // UNused
    
    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    return (
        <div>
            <div ref={containerRef} style={{ position: "relative", width: "50%", height: "100vh" }}>
                <PolygonManagerCanvas
                    handleRef={polygonManagerRef}
                    size={size}
                    transform={transform}
                    onTransform={setTransform}
                />
                <OverlyingCanvasV2
                    size={size}
                    mode="brush"
                    transform={transform}
                    onTransform={setTransform}
                    onBrush={(x, y) => polygonManagerRef.current?.addBrushAt(x, y)}
                />
            </div>
            <div style={{ width: "50%", height: "100vh", backgroundColor: "#f0f0f0" }}>

            </div>
        </div>
    )
}