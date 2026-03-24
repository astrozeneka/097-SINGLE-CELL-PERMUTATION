import { useEffect, useRef } from "react";
import { Transform } from "../utils";


export type CanvasMode = "pan" | "select" | "brush";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    mode: CanvasMode;
    transform: Transform;
    onTransform: (t: Transform) => void;
    onBrush: (x: number, y: number) => void;
    onBrushMove?: (x: number, y: number) => void;
    onCursorMove?: (x: number, y: number) => void;
    onCursorLeave?: () => void;
    onBrushResize?: (delta: number) => void;
}

export function OverlyingCanvasV2({ size, mode, transform, onTransform, onBrush, onBrushMove, onCursorMove, onCursorLeave, onBrushResize }: OverlyingCanvasParams) {
    const ref = useRef<HTMLCanvasElement>(null);
    const isDown = useRef(false);
    const onBrushResizeRef = useRef(onBrushResize);
    onBrushResizeRef.current = onBrushResize;

    useEffect(() => {
        const canvas = ref.current!;
        const onWheel = (e: WheelEvent) => {
            if (e.shiftKey) {
                e.preventDefault();
                onBrushResizeRef.current?.(e.deltaY);
            }
        };
        canvas.addEventListener("wheel", onWheel, { passive: false });
        return () => canvas.removeEventListener("wheel", onWheel);
    }, []);

    function canvasPos(e: React.MouseEvent) {
        const r = ref.current!.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    return (
        <canvas
            ref={ref}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: mode === "brush" ? "none" : mode === "pan" ? "grab" : "crosshair" }}
            onPointerDown={e => { isDown.current = true; if (mode === "brush") { const { x, y } = canvasPos(e); onBrush(x, y); } }}
            onPointerUp={() => { isDown.current = false; }}
            onMouseMove={e => {
                const { x, y } = canvasPos(e);
                onCursorMove?.(x, y);
                if (mode === "brush" && isDown.current) onBrushMove?.(x, y);
            }}
            onMouseLeave={() => onCursorLeave?.()}
        />
    );
}
