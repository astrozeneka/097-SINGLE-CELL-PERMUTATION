import { useRef } from "react";

const BRUSH_RADIUS = 20;

// Points must be in the same normalized space as the data passed to the WebGL scatter component.
// matrix is the full MVP matrix (column-major, WebGL convention), same as passed to Scatter3DV2.
export function OverlyingCanvasV3({ size, matrix, points, onBrush }: {
    size: { w: number; h: number };
    matrix: Float32Array;
    points: { x: number; y: number; z: number }[];
    onBrush: (indices: Set<number>) => void;
}) {
    const canvasRef   = useRef<HTMLCanvasElement>(null);
    const isDown      = useRef(false);
    const cursorPos   = useRef<{ x: number; y: number } | null>(null);
    const lastBrush   = useRef(0);
    const matrixRef   = useRef(matrix);
    matrixRef.current = matrix;
    const pointsRef   = useRef(points);
    pointsRef.current = points;
    const sizeRef     = useRef(size);
    sizeRef.current   = size;
    const onBrushRef  = useRef(onBrush);
    onBrushRef.current = onBrush;

    function canvasPos(e: React.MouseEvent) {
        const r = canvasRef.current!.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    // Project a normalized 3D point to screen coordinates using the MVP matrix.
    // Column-major: result.x = m[0]*x + m[4]*y + m[8]*z + m[12], etc.
    function toScreen(px: number, py: number, pz: number): { sx: number; sy: number } {
        const m = matrixRef.current;
        const { w, h } = sizeRef.current;
        const cx = m[0]*px + m[4]*py + m[8]*pz  + m[12];
        const cy = m[1]*px + m[5]*py + m[9]*pz  + m[13];
        const cw = m[3]*px + m[7]*py + m[11]*pz + m[15];
        return { sx: (cx / cw + 1) / 2 * w, sy: (1 - cy / cw) / 2 * h };
    }

    function selectWithin(bx: number, by: number): Set<number> {
        const r2  = BRUSH_RADIUS * BRUSH_RADIUS;
        const pts = pointsRef.current;
        const out = new Set<number>();
        for (let i = 0; i < pts.length; i++) {
            const { sx, sy } = toScreen(pts[i].x, pts[i].y, pts[i].z);
            const dx = bx - sx, dy = by - sy;
            if (dx * dx + dy * dy < r2) out.add(i);
        }
        return out;
    }

    function redraw() {
        const canvas = canvasRef.current;
        if (!canvas) return;
        const { w, h } = sizeRef.current;
        if (canvas.width !== w) canvas.width = w;
        if (canvas.height !== h) canvas.height = h;
        const ctx = canvas.getContext("2d")!;
        ctx.clearRect(0, 0, w, h);
        if (!cursorPos.current) return;
        const { x, y } = cursorPos.current;
        ctx.beginPath();
        ctx.arc(x, y, BRUSH_RADIUS, 0, 2 * Math.PI);
        ctx.strokeStyle = "rgba(255,255,255,0.75)";
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.stroke();
        ctx.setLineDash([]);
    }

    function handleBrush(x: number, y: number) {
        const now = Date.now();
        if (now - lastBrush.current < 16) return;
        lastBrush.current = now;
        onBrushRef.current(selectWithin(x, y));
    }

    return (
        <canvas
            ref={canvasRef}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: "none", width: "100%", height: "100%" }}
            onPointerDown={e => {
                isDown.current = true;
                const { x, y } = canvasPos(e);
                lastBrush.current = 0; // force immediate fire on first click
                handleBrush(x, y);
            }}
            onPointerUp={() => { isDown.current = false; }}
            onMouseMove={e => {
                const { x, y } = canvasPos(e);
                cursorPos.current = { x, y };
                redraw();
                if (isDown.current) handleBrush(x, y);
            }}
            onMouseLeave={() => { cursorPos.current = null; redraw(); }}
            onContextMenu={e => e.preventDefault()}
        />
    );
}
