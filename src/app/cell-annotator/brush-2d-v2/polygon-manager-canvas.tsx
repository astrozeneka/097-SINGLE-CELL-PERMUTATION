import { useEffect, useRef } from "react";
import polygonClipping from "polygon-clipping";

export interface PolygonManagerHandle {
    onBrushClick: (x: number, y: number) => void;
    onBrushMove: (x: number, y: number) => void;
    setCursorPos: (x: number, y: number) => void;
    clearCursor: () => void;
    adjustBrushRadius: (delta: number) => void;
}

interface PolygonManagerCanvasProps {
    size: { w: number; h: number };
    transform: any;
    onTransform: (transform: any) => void;
    handleRef?: React.MutableRefObject<PolygonManagerHandle | null>;
    subset: string;
}

export class Polygon {
    constructor(public id: string) {}
}

const BRUSH_RADIUS = 100;
const BRUSH_VERTS = 16;
const MIN_BRUSH_RADIUS = 10;
const MAX_BRUSH_RADIUS = 300;

export default function PolygonManagerCanvas({ size, transform, handleRef, subset }: PolygonManagerCanvasProps) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const polygonsBySubset = useRef<Map<string, { id: number; verts: { x: number; y: number }[] }[]>>(new Map());
    const subsetRef = useRef(subset);
    subsetRef.current = subset;
    const transformRef = useRef(transform);
    transformRef.current = transform;
    const nextId = useRef(1);
    const activePolygonId = useRef<number | null>(null);
    const lastMoveTime = useRef(0);
    const brushRadius = useRef(BRUSH_RADIUS);
    const cursorPos = useRef<{ x: number; y: number } | null>(null);

    function screenToData(sx: number, sy: number) {
        const { x, y, scale } = transformRef.current;
        return { x: (sx - x) / scale, y: -(sy - y) / scale };
    }

    function dataToScreen(dx: number, dy: number) {
        const { x, y, scale } = transformRef.current;
        return { x: dx * scale + x, y: -dy * scale + y };
    }

    useEffect(() => { activePolygonId.current = null; }, [subset]);

    function getPolygons() {
        let p = polygonsBySubset.current.get(subsetRef.current);
        if (!p) { p = []; polygonsBySubset.current.set(subsetRef.current, p); }
        return p;
    }

    function redraw() {
        const ctx = canvasRef.current!.getContext("2d")!;
        ctx.clearRect(0, 0, size.w, size.h);
        for (const poly of getPolygons()) {
            ctx.beginPath();
            poly.verts.forEach((v, i) => {
                const s = dataToScreen(v.x, v.y);
                i === 0 ? ctx.moveTo(s.x, s.y) : ctx.lineTo(s.x, s.y);
            });
            ctx.closePath();
            ctx.fillStyle = "rgba(100, 160, 255, 0.35)";
            ctx.fill();
            ctx.strokeStyle = "rgba(60, 120, 220, 0.8)";
            ctx.lineWidth = 1.5;
            ctx.stroke();
            const dcx = poly.verts.reduce((s, v) => s + v.x, 0) / poly.verts.length;
            const dcy = poly.verts.reduce((s, v) => s + v.y, 0) / poly.verts.length;
            const { x: cx, y: cy } = dataToScreen(dcx, dcy);
            ctx.fillStyle = "white";
            ctx.strokeStyle = "red";
            ctx.font = "bold 13px sans-serif";
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            ctx.fillText(String(poly.id), cx, cy);
        }
        if (cursorPos.current) {
            const { x, y } = cursorPos.current;
            ctx.beginPath();
            ctx.arc(x, y, brushRadius.current, 0, 2 * Math.PI);
            ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
            ctx.lineWidth = 1.5;
            ctx.setLineDash([4, 4]);
            ctx.stroke();
            ctx.setLineDash([]);
        }
    }

    function makeCircle(dx: number, dy: number): [number, number][] {
        const r = brushRadius.current / transformRef.current.scale;
        return Array.from({ length: BRUSH_VERTS }, (_, i) => {
            const angle = (2 * Math.PI * i) / BRUSH_VERTS;
            return [dx + r * Math.cos(angle), dy + r * Math.sin(angle)];
        });
    }

    function pointInPolygon(verts: { x: number; y: number }[], px: number, py: number): boolean {
        let inside = false;
        for (let i = 0, j = verts.length - 1; i < verts.length; j = i++) {
            const { x: xi, y: yi } = verts[i], { x: xj, y: yj } = verts[j];
            if ((yi > py) !== (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi)
                inside = !inside;
        }
        return inside;
    }

    function addBrushAt(sx: number, sy: number) {
        const { x, y } = screenToData(sx, sy);
        const id = nextId.current++;
        getPolygons().push({ id, verts: makeCircle(x, y).map(([px, py]) => ({ x: px, y: py })) });
        activePolygonId.current = id;
        redraw();
    }

    function moveBrushTo(sx: number, sy: number) {
        if (activePolygonId.current === null) return;
        const polys = getPolygons();
        const idx = polys.findIndex(p => p.id === activePolygonId.current);
        if (idx === -1) return;
        const poly = polys[idx];
        const { x, y } = screenToData(sx, sy);
        const result = polygonClipping.union([poly.verts.map(v => [v.x, v.y] as [number, number])], [makeCircle(x, y)]);
        if (result.length === 0) return;
        polys.splice(idx, 1,
            ...result.map((p, i) => ({ id: i === 0 ? poly.id : nextId.current++, verts: p[0].map(([px, py]) => ({ x: px, y: py })) }))
        );
        redraw();
    }

    function onBrushClick(sx: number, sy: number) {
        const { x, y } = screenToData(sx, sy);
        const hit = getPolygons().find(p => pointInPolygon(p.verts, x, y));
        if (hit) { activePolygonId.current = hit.id; moveBrushTo(sx, sy); }
        else addBrushAt(sx, sy);
    }

    function onBrushMove(sx: number, sy: number) {
        const now = Date.now();
        if (now - lastMoveTime.current < 16) return;
        lastMoveTime.current = now;
        moveBrushTo(sx, sy);
    }

    function setCursorPos(x: number, y: number) {
        cursorPos.current = { x, y };
        redraw();
    }

    function clearCursor() {
        cursorPos.current = null;
        redraw();
    }

    function adjustBrushRadius(delta: number) {
        brushRadius.current = Math.max(MIN_BRUSH_RADIUS, Math.min(MAX_BRUSH_RADIUS, brushRadius.current + delta * 0.1));
        redraw();
    }

    useEffect(() => {
        if (handleRef) handleRef.current = { onBrushClick, onBrushMove, setCursorPos, clearCursor, adjustBrushRadius };
        return () => { if (handleRef) handleRef.current = null; };
    });

    useEffect(() => { redraw(); }, [subset, size, transform]);

    return (
        <canvas ref={canvasRef} width={size.w} height={size.h} style={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%" }} />
    )
}
