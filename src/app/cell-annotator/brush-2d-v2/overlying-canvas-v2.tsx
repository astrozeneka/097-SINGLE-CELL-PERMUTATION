import { useRef } from "react";
import { Transform } from "../utils";


export type CanvasMode = "pan" | "select" | "brush";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    mode: CanvasMode;
    transform: Transform;
    onTransform: (t: Transform) => void;

}

export function OverlyingCanvasV2({ size, mode, transform, onTransform }: OverlyingCanvasParams) {
    const ref = useRef<HTMLCanvasElement>(null);
    
    <canvas
        ref={ref}
        width={size.w}
        height={size.h}
        style={{ position: "absolute", inset: 0, cursor: mode === "pan" ? "grab" : "crosshair" }}
    />

}