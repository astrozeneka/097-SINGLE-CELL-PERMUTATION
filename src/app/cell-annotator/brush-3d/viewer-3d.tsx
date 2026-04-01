"use client";
import { useEffect, useRef, useState } from "react";
import { mat4mul, perspectiveMat, rotateX, rotateY, Scatter3D, viewMat, DATA } from "./scatter-3d";
import { OverlyingCanvasV2 } from "../brush-2d-v2/overlying-canvas-v2";
import polygonClipping from "polygon-clipping";


const IDENTITY = new Float32Array([1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]);
const STEP = Math.PI / 36; // 5° per keypress
const BRUSH_RADIUS = 20;
const BRUSH_VERTS = 16;

// Pre-normalize DATA the same way Scatter3D does internally
const NORMALIZED_DATA = (() => {
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    for (const d of DATA) {
        if (d.x < minX) minX = d.x; if (d.x > maxX) maxX = d.x;
        if (d.y < minY) minY = d.y; if (d.y > maxY) maxY = d.y;
        if (d.z < minZ) minZ = d.z; if (d.z > maxZ) maxZ = d.z;
    }
    const span = Math.max(maxX - minX, maxY - minY, maxZ - minZ) / 1.6;
    const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2, cz = (minZ + maxZ) / 2;
    return DATA.map((d, i) => ({
        index: i,
        nx: (d.x - cx) / span,
        ny: (d.y - cy) / span,
        nz: (d.z - cz) / span,
    }));
})();


// --- Viewer3D ---
export default function Viewer3D() {
    const containerRef  = useRef<HTMLDivElement>(null);
    const drawCanvasRef = useRef<HTMLCanvasElement>(null);
    const [size, setSize]         = useState({ w: 0, h: 0 });
    const [rotation, setRotation] = useState(() => new Float32Array(IDENTITY));
    const rotationRef             = useRef(rotation);
    rotationRef.current           = rotation;
    const sizeRef                 = useRef(size);
    sizeRef.current               = size;

    const cursorPos      = useRef<{ x: number; y: number } | null>(null);
    const brushPolygon   = useRef<[number, number][][] | null>(null);
    const activePoly     = useRef<[number, number][] | null>(null);
    const lastMoveTime   = useRef(0);

    useEffect(() => {
        const el = containerRef.current!;
        const ro = new ResizeObserver(() => setSize({ w: el.clientWidth, h: el.clientHeight }));
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    useEffect(() => {
        const onKey = (e: KeyboardEvent) => {
            let delta: Float32Array | null = null;
            if (e.key === "ArrowLeft")  delta = rotateY(-STEP);
            if (e.key === "ArrowRight") delta = rotateY(+STEP);
            if (e.key === "ArrowUp")    delta = rotateX(-STEP);
            if (e.key === "ArrowDown")  delta = rotateX(+STEP);
            if (!delta) return;
            e.preventDefault();
            setRotation(mat4mul(delta, rotationRef.current));
        };
        window.addEventListener("keydown", onKey);
        return () => window.removeEventListener("keydown", onKey);
    }, []);

    const aspect = size.w > 0 ? size.w / size.h : 1;
    const P      = perspectiveMat(60 * Math.PI / 180, aspect, 0.1, 100.0);
    const V      = viewMat(3.0);
    const matrix = mat4mul(P, mat4mul(V, rotation));
    const matrixRef = useRef(matrix);
    matrixRef.current = matrix;

    function makeCircle(cx: number, cy: number): [number, number][] {
        return Array.from({ length: BRUSH_VERTS }, (_, i) => {
            const angle = (2 * Math.PI * i) / BRUSH_VERTS;
            return [cx + BRUSH_RADIUS * Math.cos(angle), cy + BRUSH_RADIUS * Math.sin(angle)];
        });
    }

    function redraw() {
        const canvas = drawCanvasRef.current;
        if (!canvas) return;
        const { w, h } = sizeRef.current;
        if (canvas.width !== w)  canvas.width  = w;
        if (canvas.height !== h) canvas.height = h;
        const ctx = canvas.getContext("2d")!;
        ctx.clearRect(0, 0, w, h);

        if (brushPolygon.current) {
            for (const ring of brushPolygon.current) {
                ctx.beginPath();
                ring.forEach(([x, y], i) => i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y));
                ctx.closePath();
                ctx.fillStyle = "rgba(100, 160, 255, 0.35)";
                ctx.fill();
                ctx.strokeStyle = "rgba(60, 120, 220, 0.8)";
                ctx.lineWidth = 1.5;
                ctx.stroke();
            }
        }

        if (cursorPos.current) {
            const { x, y } = cursorPos.current;
            ctx.beginPath();
            ctx.arc(x, y, BRUSH_RADIUS, 0, 2 * Math.PI);
            ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
            ctx.lineWidth = 1.5;
            ctx.setLineDash([4, 4]);
            ctx.stroke();
            ctx.setLineDash([]);
        }
    }

    function handleBrush(bx: number, by: number) {
        const m = matrixRef.current;
        const { w, h } = sizeRef.current;
        const circle = makeCircle(bx, by);

        if (activePoly.current) {
            const result = polygonClipping.union([activePoly.current], [circle]);
            activePoly.current = result[0]?.[0] ?? circle;
        } else {
            activePoly.current = circle;
        }
        brushPolygon.current = [activePoly.current];
        redraw();

        const selected = NORMALIZED_DATA.filter(({ nx, ny, nz }) => {
            const px = m[0]*nx + m[4]*ny + m[8]*nz  + m[12];
            const py = m[1]*nx + m[5]*ny + m[9]*nz  + m[13];
            const pw = m[3]*nx + m[7]*ny + m[11]*nz + m[15];
            const sx = (px / pw + 1) / 2 * w;
            const sy = (1 - py / pw) / 2 * h;
            const dx = bx - sx, dy = by - sy;
            return dx * dx + dy * dy < BRUSH_RADIUS * BRUSH_RADIUS;
        });
        console.log("selected", selected.map(d => d.index));
    }

    function handleBrushStart(bx: number, by: number) {
        activePoly.current = null;
        handleBrush(bx, by);
    }

    function handleBrushMove(bx: number, by: number) {
        const now = Date.now();
        if (now - lastMoveTime.current < 16) return;
        lastMoveTime.current = now;
        handleBrush(bx, by);
    }

    return (
        <div ref={containerRef} style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
            <Scatter3D matrix={matrix} size={size} />
            <canvas ref={drawCanvasRef} width={size.w} height={size.h}
                style={{ position: "absolute", inset: 0, width: "100%", height: "100%", pointerEvents: "none" }} />
            <OverlyingCanvasV2
                size={size}
                mode="brush"
                transform={{ x: 0, y: 0, scale: 1 }}
                onTransform={() => {}}
                onBrush={handleBrushStart}
                onBrushMove={handleBrushMove}
                onCursorMove={(x, y) => { cursorPos.current = { x, y }; redraw(); }}
                onCursorLeave={() => { cursorPos.current = null; redraw(); }}
            />
            <div style={{ position: "absolute", bottom: 16, left: 16, color: "#888", fontFamily: "monospace", fontSize: 11, pointerEvents: "none" }}>
                Arrow keys: rotate
            </div>
        </div>
    );
}
