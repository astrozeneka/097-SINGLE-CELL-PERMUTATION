import { useEffect, useRef } from "react";
import { mat4mul, rotateX, rotateY } from "../brush-3d/scatter-3d";

const DEFAULT_BRUSH_RADIUS = 20;
const ROTATION_SENSITIVITY = 0.005;

// matrix    : full MVP matrix for point projection (read-only)
// modelMatrix: the model portion only, used to compose rotation/zoom deltas
export function OverlyingCanvasV3({ size, matrix, modelMatrix, points, onBrush, onModelMatrixChange }: {
    size: { w: number; h: number };
    matrix: Float32Array;
    modelMatrix: Float32Array;
    points: { x: number; y: number; z: number }[];
    onBrush: (indices: Set<number>) => void;
    onModelMatrixChange: (m: Float32Array) => void;
}) {
    const canvasRef      = useRef<HTMLCanvasElement>(null);
    const isDown         = useRef(false);
    const isRotating     = useRef(false);
    const lastPos        = useRef({ x: 0, y: 0 });
    const cursorPos      = useRef<{ x: number; y: number } | null>(null);
    const lastBrush      = useRef(0);
    const brushRadius    = useRef(DEFAULT_BRUSH_RADIUS);

    const matrixRef               = useRef(matrix);       matrixRef.current               = matrix;
    const modelMatrixRef          = useRef(modelMatrix);  modelMatrixRef.current          = modelMatrix;
    const pointsRef               = useRef(points);       pointsRef.current               = points;
    const sizeRef                 = useRef(size);         sizeRef.current                 = size;
    const onBrushRef              = useRef(onBrush);      onBrushRef.current              = onBrush;
    const onModelMatrixChangeRef  = useRef(onModelMatrixChange); onModelMatrixChangeRef.current = onModelMatrixChange;

    // Stable redraw ref so the wheel useEffect can always call the latest version.
    const redrawRef = useRef<() => void>(() => {});

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
        ctx.arc(x, y, brushRadius.current, 0, 2 * Math.PI);
        ctx.strokeStyle = "rgba(255,255,255,0.75)";
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.stroke();
        ctx.setLineDash([]);
    }
    redrawRef.current = redraw;

    useEffect(() => {
        const canvas = canvasRef.current!;
        const onWheel = (e: WheelEvent) => {
            e.preventDefault();
            if (e.shiftKey) {
                brushRadius.current = Math.max(5, Math.min(150, brushRadius.current - e.deltaY * 0.1));
                redrawRef.current();
            } else {
                const factor = Math.pow(0.999, e.deltaY);
                const s = new Float32Array([factor,0,0,0, 0,factor,0,0, 0,0,factor,0, 0,0,0,1]);
                onModelMatrixChangeRef.current(mat4mul(s, modelMatrixRef.current) as Float32Array<ArrayBuffer>);
            }
        };
        canvas.addEventListener("wheel", onWheel, { passive: false });
        return () => canvas.removeEventListener("wheel", onWheel);
    }, []);

    function canvasPos(e: React.MouseEvent) {
        const r = canvasRef.current!.getBoundingClientRect();
        return { x: e.clientX - r.left, y: e.clientY - r.top };
    }

    function toScreen(px: number, py: number, pz: number): { sx: number; sy: number } {
        const m = matrixRef.current;
        const { w, h } = sizeRef.current;
        const cx = m[0]*px + m[4]*py + m[8]*pz  + m[12];
        const cy = m[1]*px + m[5]*py + m[9]*pz  + m[13];
        const cw = m[3]*px + m[7]*py + m[11]*pz + m[15];
        return { sx: (cx / cw + 1) / 2 * w, sy: (1 - cy / cw) / 2 * h };
    }

    function selectWithin(bx: number, by: number): Set<number> {
        const r2  = brushRadius.current * brushRadius.current;
        const pts = pointsRef.current;
        const out = new Set<number>();
        for (let i = 0; i < pts.length; i++) {
            const { sx, sy } = toScreen(pts[i].x, pts[i].y, pts[i].z);
            const dx = bx - sx, dy = by - sy;
            if (dx * dx + dy * dy < r2) out.add(i);
        }
        return out;
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
                if (e.button === 2) {
                    isRotating.current = true;
                    lastPos.current = canvasPos(e);
                } else {
                    isDown.current = true;
                    lastBrush.current = 0;
                    const { x, y } = canvasPos(e);
                    handleBrush(x, y);
                }
            }}
            onPointerUp={e => {
                if (e.button === 2) isRotating.current = false;
                else isDown.current = false;
            }}
            onMouseMove={e => {
                const { x, y } = canvasPos(e);
                cursorPos.current = { x, y };
                redraw();
                if (isRotating.current) {
                    const dx = x - lastPos.current.x;
                    const dy = y - lastPos.current.y;
                    lastPos.current = { x, y };
                    const delta = mat4mul(rotateY(dx * ROTATION_SENSITIVITY), rotateX(dy * ROTATION_SENSITIVITY));
                    onModelMatrixChangeRef.current(mat4mul(delta, modelMatrixRef.current) as Float32Array<ArrayBuffer>);
                } else if (isDown.current) {
                    handleBrush(x, y);
                }
            }}
            onMouseLeave={() => { cursorPos.current = null; redraw(); isRotating.current = false; isDown.current = false; }}
            onContextMenu={e => e.preventDefault()}
        />
    );
}
