import { useCallback, MouseEvent, useState } from "react";
import { Transform } from "./canvas";
import { Node, NodeDimensions } from "./node";

export function Selectable({ children, transform, setSelectedNodes, nodes, nodeDimensions, onContextMenu }: {
    children: React.ReactNode,
    transform: Transform,
    setSelectedNodes: (nodes: any) => void,
    nodes: Node[],
    nodeDimensions: Record<string, NodeDimensions>, // a more efficient way can be used
    onContextMenu?: (e: MouseEvent<HTMLDivElement>) => void
}) {
    

    const [selectionBox, setSelectionBox] = useState<{ startX: number; startY: number; endX: number; endY: number } | null>(null);
    const [isSelecting, setIsSelecting] = useState(false);


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

    const handleContextMenu = useCallback((e: MouseEvent<HTMLDivElement>) => {
        e.preventDefault();
        onContextMenu?.(e);
    }, [onContextMenu]);

    return (
        <div
            style={{ width: "100%", height: "100%", background: "#f5f5f5" }}
            onMouseDown={handleCanvasMouseDown}
            onMouseMove={handleCanvasMouseMove}
            onMouseUp={handleCanvasMouseUp}
            onContextMenu={handleContextMenu}
        >
            {children}

            {selectionBox && (() => {
                // Convert canvas coordinates back to screen coordinates
                const screenStartX = selectionBox.startX * transform.scale + transform.x;
                const screenStartY = selectionBox.startY * transform.scale + transform.y;
                const screenEndX = selectionBox.endX * transform.scale + transform.x;
                const screenEndY = selectionBox.endY * transform.scale + transform.y;

                return (
                    <div
                        style={{
                            position: "absolute",
                            left: `${Math.min(screenStartX, screenEndX)}px`,
                            top: `${Math.min(screenStartY, screenEndY)}px`,
                            width: `${Math.abs(screenEndX - screenStartX)}px`,
                            height: `${Math.abs(screenEndY - screenStartY)}px`,
                            border: "2px dashed #4A90E2",
                            background: "rgba(74, 144, 226, 0.1)",
                            pointerEvents: "none",
                            zIndex: 1000
                        }}
                    />
                );
            })()}
        </div>
    )
}