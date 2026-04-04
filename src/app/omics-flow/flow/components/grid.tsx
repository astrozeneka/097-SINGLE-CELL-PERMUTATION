import { NodeManager } from "../providers/node-manager";
import { NodeComponent } from "./node";
import { Canvas } from "./canvas";
import { Edge } from "./edge";

export function Grid() {
    const nodeManager = NodeManager.instance;
    const nodes = nodeManager.getNodes();
    const nodesMap = nodes.reduce((acc, node) => {
        acc[node.uid] = node;
        return acc;
    }, {} as Record<string, typeof nodes[0]>);

    const edges = nodes.flatMap(node =>
        node.edgesIn.map(sourceUid => ({
            source: nodesMap[sourceUid],
            destination: node
        }))
    ).filter(edge => edge.source);

    return (
        <div style={{ width: "100%", height: "100%", background: "#f5f5f5" }}>
            <Canvas>
                {edges.map(edge => (
                    <Edge
                        key={`${edge.source.uid}-${edge.destination.uid}`}
                        source={edge.source}
                        destination={edge.destination}
                    />
                ))}
                {nodes.map(node => (
                    <NodeComponent key={node.uid} node={node} />
                ))}
            </Canvas>
        </div>
    );
}