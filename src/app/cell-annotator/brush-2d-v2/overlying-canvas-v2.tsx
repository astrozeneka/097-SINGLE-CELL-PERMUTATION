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
    const isPanning = useRef(false);
    const lastPos = useRef({ x: 0, y: 0 });
    const onBrushResizeRef = useRef(onBrushResize);
    onBrushResizeRef.current = onBrushResize;
    const transformRef = useRef(transform);
    transformRef.current = transform;
    const onTransformRef = useRef(onTransform);
    onTransformRef.current = onTransform;

    useEffect(() => {
        console.log("WHEEL init")
        const canvas = ref.current!;
        const onWheel = (e: WheelEvent) => {
            e.preventDefault();
            if (e.shiftKey) {
                onBrushResizeRef.current?.(e.deltaY);
            } else {
                const r = canvas.getBoundingClientRect();
                const cx = e.clientX - r.left;
                const cy = e.clientY - r.top;
                const { x, y, scale } = transformRef.current;
                const factor = Math.pow(0.999, e.deltaY);
                onTransformRef.current({ x: cx - (cx - x) * factor, y: cy - (cy - y) * factor, scale: scale * factor });
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
            style={{ position: "absolute", inset: 0, cursor: mode === "brush" ? "none" : mode === "pan" ? "grab" : "crosshair", width: "100%", height: "100%" }}
            onPointerDown={e => {
                if (e.button === 2) {
                    isPanning.current = true;
                    lastPos.current = canvasPos(e);
                } else {
                    isDown.current = true;
                    if (mode === "brush") { const { x, y } = canvasPos(e); onBrush(x, y); }
                }
            }}
            onPointerUp={e => {
                if (e.button === 2) isPanning.current = false;
                else isDown.current = false;
            }}
            onMouseMove={e => {
                const { x, y } = canvasPos(e);
                onCursorMove?.(x, y);
                if (isPanning.current) {
                    const dx = x - lastPos.current.x;
                    const dy = y - lastPos.current.y;
                    lastPos.current = { x, y };
                    const t = transformRef.current;
                    onTransformRef.current({ ...t, x: t.x + dx, y: t.y + dy });
                } else if (mode === "brush" && isDown.current) {
                    onBrushMove?.(x, y);
                }
            }}
            onMouseLeave={() => onCursorLeave?.()}
            onContextMenu={e => e.preventDefault()}
        />
    );
}
