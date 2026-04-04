import { useRef, useState } from "react";
import { NodeManager } from "../providers/node-manager";

export interface Transform {
    x: number;
    y: number;
    scale: number;
}

export function Grid ({ initialTransform={ x: 0, y: 0, scale: 1 } }) {
    const gridRef = useRef<HTMLDivElement>(null);
    const [transform, setTransform] = useState<Transform>(initialTransform);
    const nodeManagerRef = useRef<NodeManager>(NodeManager.instance);

    return (
        <div ref={gridRef} style={{ width: "100%", height: "100%", position: "relative", overflow: "hidden", background: "red" }}>
        </div>
    );
}