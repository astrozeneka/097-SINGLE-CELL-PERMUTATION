"use client";
import { useRef } from "react";
import { Transform } from "./underlying-canvas";

interface OverlyingCanvasParams {
    size: { w: number; h: number };
    transform: Transform;
    onTransform: (t: Transform) => void;
}

export function OverlyingCanvas({ size, transform, onTransform }: OverlyingCanvasParams) {
    const ref = useRef<HTMLCanvasElement>(null);
    const drag = useRef<{ startX: number; startY: number; tx: number; ty: number } | null>(null);

    return (
        <canvas
            ref={ref}
            width={size.w}
            height={size.h}
            style={{ position: "absolute", inset: 0, cursor: "grab" }}
            onWheel={e => {
                // Zoom toward cursor: keeps the data point under the cursor fixed.
                // Derivation: new_pan = old_pan * f - cursor_offset_from_center * (f - 1)
                const rect = ref.current!.getBoundingClientRect();
                const dx = e.clientX - rect.left  - size.w / 2;
                const dy = e.clientY - rect.top   - size.h / 2;
                const f  = e.deltaY < 0 ? 1.1 : 1 / 1.1;
                onTransform({ scale: transform.scale * f, x: transform.x * f - dx * (f - 1), y: transform.y * f - dy * (f - 1) });
            }}
            onMouseDown={e => { drag.current = { startX: e.clientX, startY: e.clientY, tx: transform.x, ty: transform.y }; }}
            onMouseMove={e => {
                if (!drag.current) return;
                onTransform({ ...transform, x: drag.current.tx + e.clientX - drag.current.startX, y: drag.current.ty + e.clientY - drag.current.startY });
            }}
            onMouseUp={()    => { drag.current = null; }}
            onMouseLeave={()  => { drag.current = null; }}
        />
    );
}