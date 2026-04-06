import { useState, useCallback, useEffect, MouseEvent, Dispatch, SetStateAction } from "react";
import { NodeManager } from "../providers/node-manager";
import { Node, NodeComponent, NodeDimensions } from "./node";
import { Canvas, Transform } from "./canvas";
import { Edge } from "./edge";
import { Selectable } from "./selectable";
import { NodeContextMenu } from "./node-context-menu";
import { useSshCredentials } from "../../hook/use-ssh-credentials";

export function Grid({ selectedNodes, setSelectedNodes }: {
    selectedNodes: Node[],
    setSelectedNodes: Dispatch<SetStateAction<Node[]>>
}) {
    const nodeManager = NodeManager.instance;
    const credentials = useSshCredentials();
    const [areNodesLoading, setAreNodesLoading] = useState(true);
    const [, forceUpdate] = useState({});
    const nodes = nodeManager.getNodes();

    const nodesMap = nodes.reduce((acc, node) => {
        acc[node.uid] = node;
        return acc;
    }, {} as Record<string, typeof nodes[0]>);

    const [nodeDimensions, setNodeDimensions] = useState<Record<string, NodeDimensions>>({});
    const [contextMenu, setContextMenu] = useState<{ x: number; y: number; type: 'node' | 'canvas'; nodes: Node[] } | null>(null);
    const [transform, setTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });
    const [clipboard, setClipboard] = useState<{ nodes: Node[]; mode: 'clone' } | null>(null);

    const handleDimensionsChange = useCallback((uid: string, dimensions: NodeDimensions) => {
        setNodeDimensions(prev => ({
            ...prev,
            [uid]: dimensions
        }));
    }, []);

    const handleDrag = useCallback((uid: string, x: number, y: number) => {
        const node = nodesMap[uid];
        if (node) {
            node.x = x;
            node.y = y;
            forceUpdate({});
        }
    }, [nodesMap]);

    const handlePositionChanged = useCallback(async (uid: string, x: number, y: number) => {
        try {
            const response = await fetch(`http://192.168.64.3:3000/nodes/${uid}`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    linux_user: credentials.linux_user,
                    private_key: credentials.private_key,
                    'canvas-x': x,
                    'canvas-y': y
                })
            });

            if (!response.ok) {
                console.error(`Failed to update node position: ${response.statusText}`);
            }
        } catch (error) {
            console.error('Error updating node position:', error);
        }
    }, [credentials]);

    const handleNodeClick = useCallback((node: Node, e?: MouseEvent) => {
        if (e?.ctrlKey || e?.metaKey) {
            // Multi-select: toggle node in selection
            setSelectedNodes(prev =>
                prev.some(n => n.uid === node.uid)
                    ? prev.filter(n => n.uid !== node.uid)
                    : [...prev, node]
            );
        } else {
            // Single select
            setSelectedNodes([node]);
        }
    }, [setSelectedNodes]);
    const handleNodeContextMenu = useCallback((node: Node, e: MouseEvent<HTMLDivElement>) => {
        // If right-clicked node is not in selection, select only it
        let contextNodes: Node[];
        if (!selectedNodes.some(n => n.uid === node.uid)) {
            contextNodes = [node];
            setSelectedNodes([node]);
        } else {
            contextNodes = selectedNodes;
        }
        setContextMenu({ x: e.clientX, y: e.clientY, type: 'node', nodes: contextNodes });
    }, [selectedNodes, setSelectedNodes]);

    const handleCanvasContextMenu = useCallback((e: MouseEvent<HTMLDivElement>) => {
        // Only show canvas context menu if we have clipboard data
        if (clipboard) {
            setContextMenu({ x: e.clientX, y: e.clientY, type: 'canvas', nodes: [] });
        }
    }, [clipboard]);

    const handlePaste = useCallback(async (cursorPosition: { x: number; y: number }) => {
        if (!clipboard) return;

        // Calculate the center of the clipboard nodes
        const clipboardNodes = clipboard.nodes;
        const avgX = clipboardNodes.reduce((sum, node) => sum + node.x, 0) / clipboardNodes.length;
        const avgY = clipboardNodes.reduce((sum, node) => sum + node.y, 0) / clipboardNodes.length;

        // Convert cursor position to canvas coordinates
        const canvasX = (cursorPosition.x - transform.x) / transform.scale;
        const canvasY = (cursorPosition.y - transform.y) / transform.scale;

        // Calculate offset to center clipboard nodes around cursor
        const offsetX = canvasX - avgX;
        const offsetY = canvasY - avgY;

        // Clone each node
        for (const node of clipboardNodes) {
            const newX = node.x + offsetX;
            const newY = node.y + offsetY;

            await nodeManager.cloneNode(
                node.uid,
                newX,
                newY,
                credentials.linux_user,
                credentials.private_key
            );
        }

        // Force update to reflect new nodes
        forceUpdate({});
    }, [clipboard, transform, nodeManager, credentials]);

    const handleContextMenuAction = useCallback((action: 'copy' | 'cut' | 'delete' | 'paste', cursorPosition?: { x: number; y: number }) => {
        if (action === 'copy' && contextMenu) {
            console.log("Add things to clipboard", contextMenu.nodes)
            setClipboard({
                nodes: contextMenu.nodes,
                mode: 'clone'
            });
        } else if (action === 'paste' && clipboard && cursorPosition) {
            handlePaste(cursorPosition);
        }
        setContextMenu(null);
    }, [contextMenu, clipboard, handlePaste]);

    // Close context menu on click outside
    useEffect(() => {
        if (!contextMenu) return;

        const handleClick = () => setContextMenu(null);
        document.addEventListener('click', handleClick);
        return () => document.removeEventListener('click', handleClick);
    }, [contextMenu]);

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
        <Selectable
            transform={transform}
            setSelectedNodes={setSelectedNodes}
            nodes={nodes}
            nodeDimensions={nodeDimensions}
            onContextMenu={handleCanvasContextMenu}
        >
            <Canvas transform={transform} setTransform={setTransform}>
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
                        isSelected={selectedNodes.some(n => n.uid === node.uid)}
                        onDimensionsChange={handleDimensionsChange}
                        onDrag={handleDrag}
                        onPositionChanged={handlePositionChanged}
                        onClick={handleNodeClick}
                        onContextMenu={handleNodeContextMenu}
                    />
                ))}
            </Canvas>
            <NodeContextMenu
                contextMenu={contextMenu}
                onAction={handleContextMenuAction}
            />
        </Selectable>
    );
}