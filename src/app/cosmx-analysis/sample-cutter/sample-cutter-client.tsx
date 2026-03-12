"use client";

import { useState, Dispatch, SetStateAction } from "react";
import { DrawingCanvas } from "./drawing-canvas";
import { UnderlyingCanvas } from "./underlying-canvas";

export type Transform = { x: number; y: number; scale: number };
export type SetTransform = Dispatch<SetStateAction<Transform>>;

export function SampleCutterClient({ data }: { data: number[] }) {
    const [transform, setTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    return (
        <div style={{ position: "relative", width: "100%", height: "100vh" }}>
            <UnderlyingCanvas data={data} transform={transform} />
            <DrawingCanvas transform={transform} setTransform={setTransform} />
        </div>
    );
}
