"use client";
import { useEffect, useRef } from "react";
import { Transform } from "./underlying-canvas";

export type CanvasMode = "pan" | "select" | "polygon";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    mode: CanvasMode;
    transform: Transform;
    onTransform: (t: Transform) => void;
    // Rect corners in canvas pixel coordinates (origin = canvas top-left).
    onSelect: (rect: { x1: number; y1: number; x2: number; y2: number }) => void;
    // Polygon vertices in canvas pixel coordinates. Called live during draw and on finalize.
    // additive=true (Ctrl held) ORs with the existing mask instead of replacing it.
    onSelectPolygon: (points: { x: number; y: number }[], additive: boolean) => void;
}

// Store polygon vertices in base-clip space (normalized coords before transform) so they
// stay anchored to the correct data position even when the user scrolls during drawing.
function toBaseClip(px: number, py: number, s: { w: number; h: number }, t: Transform) {
    return {
        bx: (2 * px / s.w - 1 - 2 * t.x / s.w) / t.scale,
        by: (1 - 2 * py / s.h + 2 * t.y / s.h) / t.scale,
    };
}
function fromBaseClip(bx: number, by: number, s: { w: number; h: number }, t: Transform) {
    const cfx = bx * t.scale + 2 * t.x / s.w;
    const cfy = by * t.scale - 2 * t.y / s.h;
    return { x: (cfx + 1) * s.w / 2, y: (1 - cfy) * s.h / 2 };
}

export function OverlyingCanvas({ size, mode, transform, onTransform, onSelect, onSelectPolygon }: OverlyingCanvasParams) {
    const ref          = useRef<HTMLCanvasElement>(null);
    const drag         = useRef<{ startX: number; startY: number; tx: number; ty: number } | null>(null);
    const rectStart    = useRef<{ x: number; y: number } | null>(null);
    const polyPoints   = useRef<{ bx: number; by: number }[]>([]);
    const lastMousePos = useRef<{ x: number; y: number } | null>(null);

    // Redraw polygon preview when transform changes (e.g. scroll while drawing).
    useEffect(() => {
        if (mode !== "polygon" || polyPoints.current.length === 0) return;
        const mouse = lastMousePos.current;
        const screenPts = polyPoints.current.map(p => fromBaseClip(p.bx, p.by, size, transform));
        if (mouse) drawPolygonPreview(screenPts, mouse.x, mouse.y);
    }, [transform, size]);

    // Reset mode-specific state and clear overlay when switching modes.
    useEffect(() => {
        if (mode !== "select")  rectStart.current  = null;
        if (mode !== "polygon") polyPoints.current = [];
        ref.current!.getContext("2d")!.clearRect(0, 0, size.w, size.h);
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

    function drawPolygonPreview(points: { x: number; y: number }[], cx: number, cy: number) {
        const ctx = ref.current!.getContext("2d")!;
        ctx.clearRect(0, 0, size.w, size.h);
        if (points.length === 0) return;
        ctx.fillStyle   = "rgba(255, 255, 255, 0.08)";
        ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
        ctx.lineWidth   = 1;
        ctx.beginPath();
        ctx.moveTo(points[0].x, points[0].y);
        for (let i = 1; i < points.length; i++) ctx.lineTo(points[i].x, points[i].y);
        ctx.lineTo(cx, cy);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
    }

    return (
        <canvas
            ref={ref}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: mode === "pan" ? "grab" : "crosshair" }}
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
                } else if (mode === "polygon" && polyPoints.current.length > 0) {
                    const { x, y } = canvasPos(e);
                    lastMousePos.current = { x, y };
                    const screenPts = polyPoints.current.map(p => fromBaseClip(p.bx, p.by, size, transform));
                    drawPolygonPreview(screenPts, x, y);
                    // Live buffer update: need at least 2 placed points so cursor forms a triangle.
                    if (polyPoints.current.length >= 2)
                        onSelectPolygon([...screenPts, { x, y }], e.ctrlKey);
                }
            }}
            onMouseUp={()   => { drag.current = null; }}
            onMouseLeave={() => { drag.current = null; }}
            onClick={e => {
                if (mode === "select") {
                    const { x, y } = canvasPos(e);
                    if (!rectStart.current) {
                        rectStart.current = { x, y };
                    } else {
                        onSelect({ x1: rectStart.current.x, y1: rectStart.current.y, x2: x, y2: y });
                        rectStart.current = null;
                        ref.current!.getContext("2d")!.clearRect(0, 0, size.w, size.h);
                    }
                } else if (mode === "polygon") {
                    const { x, y } = canvasPos(e);
                    polyPoints.current = [...polyPoints.current, toBaseClip(x, y, size, transform)];
                }
            }}
            onDoubleClick={e => {
                if (mode !== "polygon") return;
                // The 2nd click of the double-click already pushed a duplicate point via onClick — remove it.
                const pts = polyPoints.current.slice(0, -1);
                polyPoints.current = [];
                ref.current!.getContext("2d")!.clearRect(0, 0, size.w, size.h);
                if (pts.length >= 3) onSelectPolygon(pts.map(p => fromBaseClip(p.bx, p.by, size, transform)), e.ctrlKey);
            }}
        />
    );
}
