"use client";
import { useEffect, useRef } from "react";
import { Transform } from "./underlying-canvas";

export type CanvasMode = "pan" | "select";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    mode: CanvasMode;
    transform: Transform;
    onTransform: (t: Transform) => void;
    // Rect corners in canvas pixel coordinates (origin = canvas top-left).
    onSelect: (rect: { x1: number; y1: number; x2: number; y2: number }) => void;
}

export function OverlyingCanvas({ size, mode, transform, onTransform, onSelect }: OverlyingCanvasParams) {
    const ref      = useRef<HTMLCanvasElement>(null);
    const drag     = useRef<{ startX: number; startY: number; tx: number; ty: number } | null>(null);
    const rectStart = useRef<{ x: number; y: number } | null>(null);

    // Clear 2D overlay and reset state when leaving select mode.
    useEffect(() => {
        if (mode !== "select") {
            rectStart.current = null;
            ref.current!.getContext("2d")!.clearRect(0, 0, size.w, size.h);
        }
    }, [mode]);

    function canvasPos(e: React.MouseEvent) {
        const r = ref.current!.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    function drawPreview(x1: number, y1: number, x2: number, y2: number) {
        const ctx = ref.current!.getContext("2d")!;
        ctx.clearRect(0, 0, size.w, size.h);
        ctx.fillStyle   = "rgba(255, 255, 255, 0.08)";
        ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
        ctx.lineWidth   = 1;
        ctx.fillRect(x1, y1, x2 - x1, y2 - y1);
        ctx.strokeRect(x1, y1, x2 - x1, y2 - y1);
    }

    return (
        <canvas
            ref={ref}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: mode === "select" ? "crosshair" : "grab" }}
            onWheel={e => {
                const r = ref.current!.getBoundingClientRect();
                const dx = e.clientX - r.left - size.w / 2;
                const dy = e.clientY - r.top  - size.h / 2;
                const f  = e.deltaY < 0 ? 1.1 : 1 / 1.1;
                onTransform({ scale: transform.scale * f, x: transform.x * f - dx * (f - 1), y: transform.y * f - dy * (f - 1) });
            }}
            onMouseDown={e => {
                if (mode !== "pan") return;
                drag.current = { startX: e.clientX, startY: e.clientY, tx: transform.x, ty: transform.y };
            }}
            onMouseMove={e => {
                if (mode === "pan" && drag.current) {
                    onTransform({ ...transform, x: drag.current.tx + e.clientX - drag.current.startX, y: drag.current.ty + e.clientY - drag.current.startY });
                } else if (mode === "select" && rectStart.current) {
                    const { x, y } = canvasPos(e);
                    drawPreview(rectStart.current.x, rectStart.current.y, x, y);
                }
            }}
            onMouseUp={()   => { drag.current = null; }}
            onMouseLeave={() => { drag.current = null; }}
            onClick={e => {
                if (mode !== "select") return;
                const { x, y } = canvasPos(e);
                if (!rectStart.current) {
                    rectStart.current = { x, y };
                } else {
                    onSelect({ x1: rectStart.current.x, y1: rectStart.current.y, x2: x, y2: y });
                    rectStart.current = null;
                    ref.current!.getContext("2d")!.clearRect(0, 0, size.w, size.h);
                }
            }}
        />
    );
}
