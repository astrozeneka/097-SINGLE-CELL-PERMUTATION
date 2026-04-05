import { useState, useCallback, useEffect } from "react";
import { NodeManager } from "../providers/node-manager";
import { Node, NodeComponent, NodeDimensions } from "./node";
import { Canvas } from "./canvas";
import { Edge } from "./edge";

export function Grid({ selectedNode, setSelectedNode }: { 
    selectedNode: Node | null, 
    setSelectedNode: (node: Node | null) => void 
}) {
    const nodeManager = NodeManager.instance;
    const [areNodesLoading, setAreNodesLoading] = useState(true);
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

    // Load the node from server effect
    useEffect(() => {
        const fetchNodes = async () => {
            setAreNodesLoading(true);
            const response = await fetch("http://192.168.64.3:3000/nodes");
            const data = await response.json();

            const nodeDataList = data.map((item: any) => ({
                info: {
                    ...item.info,
                    exports: item.info.exports || {}
                }
            }));

            nodeManager.setNodes(nodeDataList);
            forceUpdate({});
            setAreNodesLoading(false);
        };

        fetchNodes();
    }, []);

    if (areNodesLoading) {
        return (
            <div style={{ width: "100%", height: "100%", background: "#f5f5f5", display: "flex", alignItems: "center", justifyContent: "center" }}>
                Loading...
            </div>
        );
    }

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