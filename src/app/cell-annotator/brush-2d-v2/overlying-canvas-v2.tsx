import { useRef } from "react";
import { Transform } from "../utils";


export type CanvasMode = "pan" | "select" | "brush";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    mode: CanvasMode;
    transform: Transform;
    onTransform: (t: Transform) => void;
    onBrush: (x: number, y: number) => void;
    onBrushMove?: (x: number, y: number) => void;
}

export function OverlyingCanvasV2({ size, mode, transform, onTransform, onBrush, onBrushMove }: OverlyingCanvasParams) {
    const ref = useRef<HTMLCanvasElement>(null);
    const isDown = useRef(false);

    function canvasPos(e: React.MouseEvent) {
        const r = ref.current!.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    return (
        <canvas
            ref={ref}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: mode === "pan" ? "grab" : "crosshair" }}
            onClick={e => { if (mode === "brush") { const { x, y } = canvasPos(e); onBrush(x, y); } }}
            onPointerDown={() => { isDown.current = true; }}
            onPointerUp={() => { isDown.current = false; }}
            onMouseMove={e => { if (mode === "brush" && isDown.current) { const { x, y } = canvasPos(e); onBrushMove?.(x, y); } }}
        />
    );
}