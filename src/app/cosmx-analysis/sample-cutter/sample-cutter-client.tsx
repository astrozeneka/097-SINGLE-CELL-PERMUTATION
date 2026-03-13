"use client";

import { useState, Dispatch, SetStateAction } from "react";
import { DrawingCanvas } from "./drawing-canvas";
import { UnderlyingCanvas } from "./underlying-canvas";

export type Transform = { x: number; y: number; scale: number };
export type SetTransform = Dispatch<SetStateAction<Transform>>;
export type SavedPolygon = { name: string; points: { x: number; y: number }[] };
export type DataBounds = { minX: number; maxX: number; minY: number; maxY: number };

export function SampleCutterClient({ data, numClusters, savedPolygons, dataBounds }: { data: number[]; numClusters: number; savedPolygons: SavedPolygon[]; dataBounds: DataBounds }) {
    const [transform, setTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    return (
        <div style={{ position: "relative", width: "100%", height: "100vh" }}>
            <UnderlyingCanvas data={data} numClusters={numClusters} transform={transform} />
            <DrawingCanvas transform={transform} setTransform={setTransform} savedPolygons={savedPolygons} dataBounds={dataBounds} />
        </div>
    );
}
