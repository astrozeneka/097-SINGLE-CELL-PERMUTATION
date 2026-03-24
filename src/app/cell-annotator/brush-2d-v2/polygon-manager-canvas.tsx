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

export default function PolygonManagerCanvas({ size, handleRef, subset }: PolygonManagerCanvasProps) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const polygonsBySubset = useRef<Map<string, { id: number; verts: { x: number; y: number }[] }[]>>(new Map());
    const subsetRef = useRef(subset);
    subsetRef.current = subset;
    const nextId = useRef(1);
    const activePolygonId = useRef<number | null>(null);
    const lastMoveTime = useRef(0);
    const brushRadius = useRef(BRUSH_RADIUS);
    const cursorPos = useRef<{ x: number; y: number } | null>(null);

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
            poly.verts.forEach((v, i) => i === 0 ? ctx.moveTo(v.x, v.y) : ctx.lineTo(v.x, v.y));
            ctx.closePath();
            ctx.fillStyle = "rgba(100, 160, 255, 0.35)";
            ctx.fill();
            ctx.strokeStyle = "rgba(60, 120, 220, 0.8)";
            ctx.lineWidth = 1.5;
            ctx.stroke();
            const cx = poly.verts.reduce((s, v) => s + v.x, 0) / poly.verts.length;
            const cy = poly.verts.reduce((s, v) => s + v.y, 0) / poly.verts.length;
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

    function makeCircle(x: number, y: number): [number, number][] {
        return Array.from({ length: BRUSH_VERTS }, (_, i) => {
            const angle = (2 * Math.PI * i) / BRUSH_VERTS;
            return [x + brushRadius.current * Math.cos(angle), y + brushRadius.current * Math.sin(angle)];
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

    function addBrushAt(x: number, y: number) {
        const id = nextId.current++;
        getPolygons().push({ id, verts: makeCircle(x, y).map(([px, py]) => ({ x: px, y: py })) });
        activePolygonId.current = id;
        redraw();
    }

    function moveBrushTo(x: number, y: number) {
        if (activePolygonId.current === null) return;
        const polys = getPolygons();
        const idx = polys.findIndex(p => p.id === activePolygonId.current);
        if (idx === -1) return;
        const poly = polys[idx];
        const result = polygonClipping.union([poly.verts.map(v => [v.x, v.y] as [number, number])], [makeCircle(x, y)]);
        if (result.length === 0) return;
        polys.splice(idx, 1,
            ...result.map((p, i) => ({ id: i === 0 ? poly.id : nextId.current++, verts: p[0].map(([px, py]) => ({ x: px, y: py })) }))
        );
        redraw();
    }

    function onBrushClick(x: number, y: number) {
        const hit = getPolygons().find(p => pointInPolygon(p.verts, x, y));
        if (hit) { activePolygonId.current = hit.id; moveBrushTo(x, y); }
        else addBrushAt(x, y);
    }

    function onBrushMove(x: number, y: number) {
        console.log("Brush move to", x, y);
        const now = Date.now();
        if (now - lastMoveTime.current < 16) return;
        lastMoveTime.current = now;
        moveBrushTo(x, y);
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

    useEffect(() => { redraw(); }, [subset, size]);

    return (
        <canvas ref={canvasRef} width={size.w} height={size.h} style={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%" }} />
    )
}
