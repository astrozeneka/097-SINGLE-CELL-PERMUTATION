"use client";

import { Children, useRef, useState } from "react";

interface HorizontalSplitProps {
    children: React.ReactNode;
    initialPanelWidth?: number;
    minPanelWidth?: number;
    maxPanelWidth?: number;
}

export default function HorizontalSplit({
    children,
    initialPanelWidth = 500,
    minPanelWidth = 200,
    maxPanelWidth = 900,
}: HorizontalSplitProps) {
    const panels = Children.toArray(children);
    const [widths, setWidths] = useState<number[]>(() => panels.slice(1).map(() => initialPanelWidth));
    const dragRef = useRef<{ startX: number; startW: number; index: number } | null>(null);

    const onDividerMouseDown = (e: React.MouseEvent, index: number) => {
        e.preventDefault();
        dragRef.current = { startX: e.clientX, startW: widths[index], index };
        const onMove = (ev: MouseEvent) => {
            const { startX, startW, index } = dragRef.current!;
            const delta = startX - ev.clientX;
            const newW = Math.max(minPanelWidth, Math.min(maxPanelWidth, startW + delta));
            setWidths(ws => ws.map((w, i) => i === index ? newW : w));
        };
        const onUp = () => {
            dragRef.current = null;
            window.removeEventListener("mousemove", onMove);
            window.removeEventListener("mouseup", onUp);
        };
        window.addEventListener("mousemove", onMove);
        window.addEventListener("mouseup", onUp);
    };

    return (
        <div style={{ display: "flex", flexDirection: "row", width: "100%", height: "100%" }}>
            {panels.map((panel, i) => (
                <>
                    {i > 0 && (
                        <div
                            key={`divider-${i}`}
                            onMouseDown={e => onDividerMouseDown(e, i - 1)}
                            style={{ width: 5, flexShrink: 0, cursor: "col-resize", background: "#c8c8c8", transition: "background 0.15s ease" }}
                            onMouseEnter={e => (e.currentTarget.style.background = "#0078d4")}
                            onMouseLeave={e => (e.currentTarget.style.background = "#c8c8c8")}
                        />
                    )}
                    <div
                        key={`panel-${i}`}
                        style={i === 0
                            ? { flex: "1 1 auto", position: "relative", overflow: "hidden" }
                            : { flex: `0 0 ${widths[i - 1]}px`, position: "relative", overflow: "hidden" }
                        }
                    >
                        {panel}
                    </div>
                </>
            ))}
        </div>
    );
}
