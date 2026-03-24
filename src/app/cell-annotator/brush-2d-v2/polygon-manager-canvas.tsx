import { useEffect, useRef } from "react";
import polygonClipping from "polygon-clipping";

export interface PolygonManagerHandle {
    onBrushClick: (x: number, y: number) => void;
    onBrushMove: (x: number, y: number) => void;
}

interface PolygonManagerCanvasProps {
    size: { w: number; h: number };
    transform: any;
    onTransform: (transform: any) => void;
    handleRef?: React.MutableRefObject<PolygonManagerHandle | null>;
}

export class Polygon {
    constructor(public id: string) {}
}

const BRUSH_RADIUS = 100;
const BRUSH_VERTS = 16;

export default function PolygonManagerCanvas({ size, handleRef }: PolygonManagerCanvasProps) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const polygons = useRef<{ id: number; verts: { x: number; y: number }[] }[]>([]);
    const nextId = useRef(1);
    const activePolygonId = useRef<number | null>(null);
    const lastMoveTime = useRef(0);

    function redraw() {
        const ctx = canvasRef.current!.getContext("2d")!;
        ctx.clearRect(0, 0, size.w, size.h);
        for (const poly of polygons.current) {
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
        console.log("==", polygons.current)
    }

    function makeCircle(x: number, y: number): [number, number][] {
        return Array.from({ length: BRUSH_VERTS }, (_, i) => {
            const angle = (2 * Math.PI * i) / BRUSH_VERTS;
            return [x + BRUSH_RADIUS * Math.cos(angle), y + BRUSH_RADIUS * Math.sin(angle)];
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
        polygons.current.push({ id, verts: makeCircle(x, y).map(([px, py]) => ({ x: px, y: py })) });
        activePolygonId.current = id;
        redraw();
        console.log("Added polygon at", { x, y });
    }

    function moveBrushTo(x: number, y: number) {
        if (activePolygonId.current === null) return;
        const idx = polygons.current.findIndex(p => p.id === activePolygonId.current);
        if (idx === -1) return;
        const poly = polygons.current[idx];
        const result = polygonClipping.union([poly.verts.map(v => [v.x, v.y] as [number, number])], [makeCircle(x, y)]);
        polygons.current.splice(idx, 1,
            ...result.map((p, i) => ({ id: i === 0 ? poly.id : nextId.current++, verts: p[0].map(([px, py]) => ({ x: px, y: py })) }))
        );
        redraw();
    }

    function onBrushClick(x: number, y: number) {
        const hit = polygons.current.find(p => pointInPolygon(p.verts, x, y));
        if (hit) { activePolygonId.current = hit.id; moveBrushTo(x, y); }
        else addBrushAt(x, y);
    }

    function onBrushMove(x: number, y: number) {
        const now = Date.now();
        if (now - lastMoveTime.current < 16) return;
        lastMoveTime.current = now;
        moveBrushTo(x, y);
    }

    useEffect(() => {
        if (handleRef) handleRef.current = { onBrushClick, onBrushMove };
        return () => { if (handleRef) handleRef.current = null; };
    });

    return (
        <canvas ref={canvasRef} width={size.w} height={size.h} style={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%" }} />
    )
}
