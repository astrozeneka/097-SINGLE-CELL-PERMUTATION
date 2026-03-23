import { useEffect, useRef } from "react";

export interface PolygonManagerHandle {
    addBrushAt: (x: number, y: number) => void;
}

interface PolygonManagerCanvasProps {
    size: { w: number; h: number };
    transform: any;
    onTransform: (transform: any) => void;
    handleRef?: React.MutableRefObject<PolygonManagerHandle | null>;
}

const BRUSH_RADIUS = 100;
const BRUSH_VERTS = 16;

export default function PolygonManagerCanvas({ size, handleRef }: PolygonManagerCanvasProps) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const polygons = useRef<{ id: number; verts: { x: number; y: number }[] }[]>([]);
    const nextId = useRef(1);

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

    function addBrushAt(x: number, y: number) {
        const verts = Array.from({ length: BRUSH_VERTS }, (_, i) => {
            const angle = (2 * Math.PI * i) / BRUSH_VERTS;
            return { x: x + BRUSH_RADIUS * Math.cos(angle), y: y + BRUSH_RADIUS * Math.sin(angle) };
        });
        polygons.current.push({ id: nextId.current++, verts });
        redraw();
        console.log("Added polygon at", { x, y }, "with verts", verts);
    }

    useEffect(() => {
        if (handleRef) handleRef.current = { addBrushAt };
        return () => { if (handleRef) handleRef.current = null; };
    });

    return (
        <canvas ref={canvasRef} width={size.w} height={size.h} style={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%" }} />
    )
}
