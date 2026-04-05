import { useState, useCallback, useEffect, MouseEvent, Dispatch, SetStateAction } from "react";
import { NodeManager } from "../providers/node-manager";
import { Node, NodeComponent, NodeDimensions } from "./node";
import { Canvas, Transform } from "./canvas";
import { Edge } from "./edge";

export function Grid({ selectedNodes, setSelectedNodes }: {
    selectedNodes: Node[],
    setSelectedNodes: Dispatch<SetStateAction<Node[]>>
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
    const [selectionBox, setSelectionBox] = useState<{ startX: number; startY: number; endX: number; endY: number } | null>(null);
    const [isSelecting, setIsSelecting] = useState(false);
    const [contextMenu, setContextMenu] = useState<{ x: number; y: number } | null>(null);
    const [transform, setTransform] = useState<Transform>({ x: 0, y: 0, scale: 1 });

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

    const handleCanvasMouseDown = useCallback((e: MouseEvent<HTMLDivElement>) => {
        // Only start selection on left mouse button
        // Node clicks will stopPropagation, so this only fires when clicking canvas background
        if (e.button === 0) {
            const rect = (e.currentTarget as HTMLElement).getBoundingClientRect();
            const screenX = e.clientX - rect.left;
            const screenY = e.clientY - rect.top;

            // Convert screen coordinates to canvas coordinates
            const canvasX = (screenX - transform.x) / transform.scale;
            const canvasY = (screenY - transform.y) / transform.scale;

            setSelectionBox({ startX: canvasX, startY: canvasY, endX: canvasX, endY: canvasY });
            setIsSelecting(true);
        }
    }, [transform]);

    const handleCanvasMouseMove = useCallback((e: MouseEvent<HTMLDivElement>) => {
        if (isSelecting && selectionBox) {
            const rect = (e.currentTarget as HTMLElement).getBoundingClientRect();
            const screenX = e.clientX - rect.left;
            const screenY = e.clientY - rect.top;

            // Convert screen coordinates to canvas coordinates
            const canvasX = (screenX - transform.x) / transform.scale;
            const canvasY = (screenY - transform.y) / transform.scale;

            setSelectionBox({ ...selectionBox, endX: canvasX, endY: canvasY });
        }
    }, [isSelecting, selectionBox, transform]);

    const handleCanvasMouseUp = useCallback(() => {
        if (isSelecting && selectionBox) {
            // Calculate selection box bounds
            const minX = Math.min(selectionBox.startX, selectionBox.endX);
            const maxX = Math.max(selectionBox.startX, selectionBox.endX);
            const minY = Math.min(selectionBox.startY, selectionBox.endY);
            const maxY = Math.max(selectionBox.startY, selectionBox.endY);

            // Find nodes that intersect with selection box
            const selectedNodesList = nodes.filter(node => {
                const dims = nodeDimensions[node.uid];
                if (!dims) return false;

                // Check if any corner of the node is inside the selection box
                const nodeMinX = node.x;
                const nodeMaxX = node.x + dims.width;
                const nodeMinY = node.y;
                const nodeMaxY = node.y + dims.height;

                // Any corner inside selection box
                const topLeftInside = nodeMinX >= minX && nodeMinX <= maxX && nodeMinY >= minY && nodeMinY <= maxY;
                const topRightInside = nodeMaxX >= minX && nodeMaxX <= maxX && nodeMinY >= minY && nodeMinY <= maxY;
                const bottomLeftInside = nodeMinX >= minX && nodeMinX <= maxX && nodeMaxY >= minY && nodeMaxY <= maxY;
                const bottomRightInside = nodeMaxX >= minX && nodeMaxX <= maxX && nodeMaxY >= minY && nodeMaxY <= maxY;

                // For now, use any corner inside logic
                return topLeftInside || topRightInside || bottomLeftInside || bottomRightInside;

                // Alternative: All corners inside (commented out as per requirement)
                // const allCornersInside =
                //     (nodeMinX >= minX && nodeMinX <= maxX && nodeMinY >= minY && nodeMinY <= maxY) &&
                //     (nodeMaxX >= minX && nodeMaxX <= maxX && nodeMinY >= minY && nodeMinY <= maxY) &&
                //     (nodeMinX >= minX && nodeMinX <= maxX && nodeMaxY >= minY && nodeMaxY <= maxY) &&
                //     (nodeMaxX >= minX && nodeMaxX <= maxX && nodeMaxY >= minY && nodeMaxY <= maxY);
                // return allCornersInside;
            });

            setSelectedNodes(selectedNodesList);
            setSelectionBox(null);
            setIsSelecting(false);
        }
    }, [isSelecting, selectionBox, nodes, nodeDimensions, setSelectedNodes]);

    const handleNodeContextMenu = useCallback((node: Node, e: MouseEvent<HTMLDivElement>) => {
        // If right-clicked node is not in selection, select only it
        if (!selectedNodes.some(n => n.uid === node.uid)) {
            setSelectedNodes([node]);
        }
        setContextMenu({ x: e.clientX, y: e.clientY });
    }, [selectedNodes, setSelectedNodes]);

    const handleContextMenuAction = useCallback((action: 'copy' | 'cut' | 'delete') => {
        // Placeholder implementations
        console.log(`${action} action on nodes:`, selectedNodes.map(n => n.uid));
        setContextMenu(null);
    }, [selectedNodes]);

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
        <div
            style={{ width: "100%", height: "100%", background: "#f5f5f5" }}
            onMouseDown={handleCanvasMouseDown}
            onMouseMove={handleCanvasMouseMove}
            onMouseUp={handleCanvasMouseUp}
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
                        onPositionChange={handlePositionChange}
                        onClick={handleNodeClick}
                        onContextMenu={handleNodeContextMenu}
                    />
                ))}
                {selectionBox && (
                    <div
                        style={{
                            position: "absolute",
                            left: `${Math.min(selectionBox.startX, selectionBox.endX)}px`,
                            top: `${Math.min(selectionBox.startY, selectionBox.endY)}px`,
                            width: `${Math.abs(selectionBox.endX - selectionBox.startX)}px`,
                            height: `${Math.abs(selectionBox.endY - selectionBox.startY)}px`,
                            border: "2px dashed #4A90E2",
                            background: "rgba(74, 144, 226, 0.1)",
                            pointerEvents: "none",
                            zIndex: 1000
                        }}
                    />
                )}
            </Canvas>
            {contextMenu && (
                <div
                    style={{
                        position: "fixed",
                        left: `${contextMenu.x}px`,
                        top: `${contextMenu.y}px`,
                        background: "white",
                        border: "1px solid #ccc",
                        borderRadius: "4px",
                        boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
                        zIndex: 10000,
                        minWidth: "150px"
                    }}
                >
                    <div
                        onClick={() => handleContextMenuAction('copy')}
                        style={{
                            padding: "8px 16px",
                            cursor: "pointer",
                            borderBottom: "1px solid #eee"
                        }}
                        onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                        onMouseLeave={(e) => e.currentTarget.style.background = "white"}
                    >
                        Copy clone
                    </div>
                    <div
                        onClick={() => handleContextMenuAction('cut')}
                        style={{
                            padding: "8px 16px",
                            cursor: "pointer",
                            borderBottom: "1px solid #eee"
                        }}
                        onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                        onMouseLeave={(e) => e.currentTarget.style.background = "white"}
                    >
                        Cut
                    </div>
                    <div
                        onClick={() => handleContextMenuAction('delete')}
                        style={{
                            padding: "8px 16px",
                            cursor: "pointer",
                            color: "#d32f2f"
                        }}
                        onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                        onMouseLeave={(e) => e.currentTarget.style.background = "white"}
                    >
                        Delete
                    </div>
                </div>
            )}
        </div>
    );
}