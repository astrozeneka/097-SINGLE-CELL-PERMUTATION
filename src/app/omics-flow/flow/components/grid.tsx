import { NodeManager } from "../providers/node-manager";
import { NodeComponent } from "./node";
import { Canvas } from "./canvas";

export function Grid() {
    const nodeManager = NodeManager.instance;
    const nodes = nodeManager.getNodes();

    return (
        <div style={{ width: "100%", height: "100%", background: "#f5f5f5" }}>
            <Canvas>
                {nodes.map(node => (
                    <NodeComponent key={node.uid} node={node} />
                ))}
            </Canvas>
        </div>
    );
}