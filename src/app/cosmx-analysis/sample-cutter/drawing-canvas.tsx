"use client";

import { useRef, useEffect, MouseEvent } from "react";
import { Transform, SetTransform } from "./sample-cutter-client";

export function DrawingCanvas({ transform, setTransform }: { transform: Transform; setTransform: SetTransform }) {
    const canvasRef = useRef<HTMLCanvasElement>(null);
    const dragRef = useRef<{ startX: number; startY: number; origX: number; origY: number } | null>(null);

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
        return () => canvas.removeEventListener("wheel", onWheel);
    }, [setTransform]);

    useEffect(() => {
        const canvas = canvasRef.current!;
        const ctx = canvas.getContext("2d")!;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.save();
        ctx.translate(transform.x, transform.y);
        ctx.scale(transform.scale, transform.scale);
        ctx.strokeStyle = "red";
        ctx.lineWidth = 2 / transform.scale;
        ctx.strokeRect(100, 100, 200, 150);
        ctx.restore();
    }, [transform]);

    function onMouseDown(e: MouseEvent) {
        dragRef.current = { startX: e.clientX, startY: e.clientY, origX: transform.x, origY: transform.y };
    }

    function onMouseMove(e: MouseEvent) {
        if (!dragRef.current) return;
        const { startX, startY, origX, origY } = dragRef.current;
        const dx = e.clientX - startX;
        const dy = e.clientY - startY;
        setTransform(t => ({ ...t, x: origX + dx, y: origY + dy }));
    }

    function onMouseUp() { dragRef.current = null; }

    function onClick(_e: MouseEvent) {
        // TODO: handle polygon point placement
    }

    return (
        <canvas
            ref={canvasRef}
            style={{ position: "absolute", inset: 0, width: "100%", height: "100%", cursor: "crosshair" }}
            onMouseDown={onMouseDown}
            onMouseMove={onMouseMove}
            onMouseUp={onMouseUp}
            onMouseLeave={onMouseUp}
            onClick={onClick}
        />
    );
}
