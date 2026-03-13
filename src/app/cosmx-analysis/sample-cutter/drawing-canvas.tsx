"use client";

import { useRef, useEffect, useState, MouseEvent } from "react";
import { Transform, SetTransform, SavedPolygon, DataBounds } from "./sample-cutter-client";

type Point = { x: number; y: number };

export function DrawingCanvas({ transform, setTransform, savedPolygons, dataBounds }: { transform: Transform; setTransform: SetTransform; savedPolygons: SavedPolygon[]; dataBounds: DataBounds }) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const dragRef = useRef<{ startX: number; startY: number; origX: number; origY: number } | null>(null);
    const [currentPoly, setCurrentPoly] = useState<Point[]>([]);
    const [polygons, setPolygons] = useState<Point[][]>([]);
    const [copied, setCopied] = useState(false);
    const currentPolyRef = useRef(currentPoly);
    currentPolyRef.current = currentPoly;

    useEffect(() => {
        const canvas = canvasRef.current!;
        canvas.width = canvas.clientWidth;
        canvas.height = canvas.clientHeight;

        const onWheel = (e: WheelEvent) => {
            e.preventDefault();
            const factor = e.deltaY < 0 ? 1.1 : 1 / 1.1;
            const rect = canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            setTransform(t => ({
                scale: t.scale * factor,
                x: mx - (mx - t.x) * factor,
                y: my - (my - t.y) * factor,
            }));
        };
        canvas.addEventListener("wheel", onWheel, { passive: false });

        const onKeyDown = (e: KeyboardEvent) => {
            if (e.key !== "c" && e.key !== "C") return;
            const pts = currentPolyRef.current;
            if (pts.length < 2) return;
            setPolygons(ps => [...ps, pts]);
            setCurrentPoly([]);
            const ring = [...pts, pts[0]].map(p => [p.x, p.y]);
            navigator.clipboard.writeText(JSON.stringify({ type: "Feature", geometry: { type: "Polygon", coordinates: [ring] } }));
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
        };
        window.addEventListener("keydown", onKeyDown);

        return () => {
            canvas.removeEventListener("wheel", onWheel);
            window.removeEventListener("keydown", onKeyDown);
        };
    }, [setTransform]);

    useEffect(() => {
        const canvas = canvasRef.current!;
        const ctx = canvas.getContext("2d")!;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.save();
        ctx.translate(transform.x, transform.y);
        ctx.scale(transform.scale, transform.scale);
        ctx.lineWidth = 2 / transform.scale;

        const { minX, maxX, minY, maxY } = dataBounds;
        const midX = (minX + maxX) / 2;
        const midY = (minY + maxY) / 2;
        const W = canvas.width, H = canvas.height;
        const s = Math.min(W, H) / (maxX - minX);
        const toPixel = (p: Point) => ({
            x: W / 2 + (p.x - midX) * s,
            y: H / 2 - (p.y - midY) * s,
        });

        const stroke = (pts: Point[], close: boolean) => {
            if (pts.length < 2) return;
            const mapped = pts.map(toPixel);
            ctx.beginPath();
            ctx.moveTo(mapped[0].x, mapped[0].y);
            for (let i = 1; i < mapped.length; i++) ctx.lineTo(mapped[i].x, mapped[i].y);
            if (close) ctx.closePath();
            ctx.stroke();
        };

        ctx.strokeStyle = "cyan";
        ctx.fillStyle = "cyan";
        ctx.font = `${12 / transform.scale}px sans-serif`;
        for (const { name, points } of savedPolygons) {
            stroke(points, true);
            if (points.length > 0) { const p = toPixel(points[0]); ctx.fillText(name, p.x, p.y); }
        }

        ctx.strokeStyle = "yellow";
        for (const poly of polygons) stroke(poly, true);

        ctx.strokeStyle = "red";
        stroke(currentPoly, false);

        ctx.restore();
    }, [transform, currentPoly, polygons, savedPolygons, dataBounds]);

    function onMouseDown(e: MouseEvent) {
        if (e.button !== 2) return;
        dragRef.current = { startX: e.clientX, startY: e.clientY, origX: transform.x, origY: transform.y };
    }

    function onMouseMove(e: MouseEvent) {
        if (!dragRef.current) return;
        const { startX, startY, origX, origY } = dragRef.current;
        setTransform(t => ({ ...t, x: origX + (e.clientX - startX), y: origY + (e.clientY - startY) }));
    }

    function onMouseUp(e: MouseEvent) {
        if (e.button === 2) dragRef.current = null;
    }

    function onClick(e: MouseEvent) {
        const canvas = canvasRef.current!;
        const rect = canvas.getBoundingClientRect();
        const baseX = (e.clientX - rect.left - transform.x) / transform.scale;
        const baseY = (e.clientY - rect.top - transform.y) / transform.scale;
        const { minX, maxX, minY, maxY } = dataBounds;
        const midX = (minX + maxX) / 2;
        const midY = (minY + maxY) / 2;
        const s = Math.min(canvas.width, canvas.height) / (maxX - minX);
        const x = midX + (baseX - canvas.width / 2) / s;
        const y = midY - (baseY - canvas.height / 2) / s;
        setCurrentPoly(pts => [...pts, { x, y }]);
    }

    return (
        <>
        {copied && (
            <div style={{ position: "absolute", bottom: 16, left: "50%", transform: "translateX(-50%)", background: "rgba(0,0,0,0.75)", color: "#fff", padding: "4px 14px", borderRadius: 4, fontSize: 12, pointerEvents: "none", zIndex: 10 }}>
                GeoJSON copied
            </div>
        )}
        <canvas
            ref={canvasRef}
            style={{ position: "absolute", inset: 0, width: "100%", height: "100%", cursor: "crosshair" }}
            onMouseDown={onMouseDown}
            onMouseMove={onMouseMove}
            onMouseUp={onMouseUp}
            onMouseLeave={() => { dragRef.current = null; }}
            onClick={onClick}
            onContextMenu={e => e.preventDefault()}
        />
        </>
    );
}
