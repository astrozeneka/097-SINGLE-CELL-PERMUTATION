import { useRef, useState } from "react";
import { NodeManager } from "../providers/node-manager";
import { NodeComponent } from "./node";

export interface Transform {
    x: number;
    y: number;
    scale: number;
}

export function Grid ({ initialTransform={ x: 0, y: 0, scale: 1 } }) {
    const gridRef = useRef<HTMLDivElement>(null);
    const [transform, setTransform] = useState<Transform>(initialTransform);
    const nodeManager = NodeManager.instance;

    const nodes = nodeManager.getNodes();

    return (
        <div ref={gridRef} style={{ width: "100%", height: "100%", position: "relative", overflow: "hidden", background: "#f5f5f5" }}>
            <div
                style={{
                    transform: `translate(${transform.x}px, ${transform.y}px) scale(${transform.scale})`,
                    transformOrigin: "0 0",
                    width: "100%",
                    height: "100%",
                    position: "relative",
                }}
            >
                {nodes.map(node => (
                    <NodeComponent key={node.uid} node={node} />
                ))}
            </div>
        </div>
    );
}