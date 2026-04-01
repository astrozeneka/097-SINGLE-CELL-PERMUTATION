"use client";
import { useEffect, useRef, useState } from "react";
import { mat4mul, perspectiveMat, rotateX, rotateY, Scatter3D, viewMat } from "./scatter-3d";


const IDENTITY = new Float32Array([1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]);
const STEP = Math.PI / 36; // 5° per keypress


// --- Viewer3D ---
export default function Viewer3D() {
    const containerRef = useRef<HTMLDivElement>(null);
    const [size, setSize]         = useState({ w: 0, h: 0 });
    const [rotation, setRotation] = useState(() => new Float32Array(IDENTITY));
    const rotationRef             = useRef(rotation);
    rotationRef.current           = rotation;

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

    return (
        <div ref={containerRef} style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
            <Scatter3D matrix={matrix} size={size} />
            <div style={{ position: "absolute", bottom: 16, left: 16, color: "#888", fontFamily: "monospace", fontSize: 11, pointerEvents: "none" }}>
                Arrow keys: rotate
            </div>
        </div>
    );
}
