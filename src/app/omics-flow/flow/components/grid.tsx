import { useState, useCallback } from "react";
import { NodeManager } from "../providers/node-manager";
import { NodeComponent, NodeDimensions } from "./node";
import { Canvas } from "./canvas";
import { Edge } from "./edge";

export function Grid() {
    const nodeManager = NodeManager.instance;
    const [, forceUpdate] = useState({});
    const nodes = nodeManager.getNodes();
    const nodesMap = nodes.reduce((acc, node) => {
        acc[node.uid] = node;
        return acc;
    }, {} as Record<string, typeof nodes[0]>);

    const [nodeDimensions, setNodeDimensions] = useState<Record<string, NodeDimensions>>({});

    const handleDimensionsChange = useCallback((uid: string, dimensions: NodeDimensions) => {
        setNodeDimensions(prev => ({
            ...prev,
            [uid]: dimensions
        }));
    }, []);

    const handlePositionChange = useCallback((uid: string, x: number, y: number) => {
        const node = nodesMap[uid];
        if (node) {
            node.x = x;
            node.y = y;
            forceUpdate({});
        }
    }, [nodesMap]);

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
                        sourceDimensions={nodeDimensions[edge.source.uid]}
                        destinationDimensions={nodeDimensions[edge.destination.uid]}
                    />
                ))}
                {nodes.map(node => (
                    <NodeComponent
                        key={node.uid}
                        node={node}
                        onDimensionsChange={handleDimensionsChange}
                        onPositionChange={handlePositionChange}
                    />
                ))}
            </Canvas>
        </div>
    );
}