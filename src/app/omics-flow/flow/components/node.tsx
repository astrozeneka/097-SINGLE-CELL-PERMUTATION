export interface NodeInfo {
    uid: string;
    name: string;
    description: string;
    "canvas-x": number;
    "canvas-y": number;
    "edges-in": string[];
    exports: Record<string, string>;
}

export interface NodeData {
    info: NodeInfo;
}

import { useRef, useEffect, useState, useCallback, MouseEvent } from "react";

export class Node {
    uid: string;
    name: string;
    description: string;
    x: number;
    y: number;
    edgesIn: string[];
    exports: Record<string, string>;

    constructor(data: NodeData) {
        this.uid = data.info.uid;
        this.name = data.info.name;
        this.description = data.info.description;
        this.x = data.info["canvas-x"];
        this.y = data.info["canvas-y"];
        this.edgesIn = data.info["edges-in"];
        this.exports = data.info.exports;
    }
}

export interface NodeDimensions {
    width: number;
    height: number;
}

interface NodeComponentProps {
    node: Node;
    onDimensionsChange?: (uid: string, dimensions: NodeDimensions) => void;
    onPositionChange?: (uid: string, x: number, y: number) => void;
}

type ResizeCorner = "top-left" | "top-right" | "bottom-left" | "bottom-right" | null;

export function NodeComponent({ node, onDimensionsChange, onPositionChange }: NodeComponentProps) {
    const nodeRef = useRef<HTMLDivElement>(null);
    const [isDragging, setIsDragging] = useState(false);
    const [isResizing, setIsResizing] = useState<ResizeCorner>(null);
    const [dragStart, setDragStart] = useState({ x: 0, y: 0 });
    const [resizeStart, setResizeStart] = useState({ x: 0, y: 0, width: 0, height: 0, nodeX: 0, nodeY: 0 });

    useEffect(() => {
        const element = nodeRef.current;
        if (!element || !onDimensionsChange) return;

        const updateDimensions = () => {
            const rect = element.getBoundingClientRect();
            onDimensionsChange(node.uid, {
                width: rect.width,
                height: rect.height
            });
        };

        updateDimensions();

        const resizeObserver = new ResizeObserver(updateDimensions);
        resizeObserver.observe(element);

        return () => resizeObserver.disconnect();
    }, [node.uid, onDimensionsChange]);

    useEffect(() => {
        if (!isDragging && !isResizing) return;

        const handleMouseMove = (e: globalThis.MouseEvent) => {
            if (isDragging && onPositionChange) {
                const newX = e.clientX - dragStart.x;
                const newY = e.clientY - dragStart.y;
                onPositionChange(node.uid, newX, newY);
            } else if (isResizing && nodeRef.current) {
                const element = nodeRef.current;
                const deltaX = e.clientX - resizeStart.x;
                const deltaY = e.clientY - resizeStart.y;

                let newWidth = resizeStart.width;
                let newHeight = resizeStart.height;
                let newX = resizeStart.nodeX;
                let newY = resizeStart.nodeY;

                switch (isResizing) {
                    case "top-left":
                        newWidth = resizeStart.width - deltaX;
                        newHeight = resizeStart.height - deltaY;
                        newX = resizeStart.nodeX + deltaX;
                        newY = resizeStart.nodeY + deltaY;
                        break;
                    case "top-right":
                        newWidth = resizeStart.width + deltaX;
                        newHeight = resizeStart.height - deltaY;
                        newY = resizeStart.nodeY + deltaY;
                        break;
                    case "bottom-left":
                        newWidth = resizeStart.width - deltaX;
                        newHeight = resizeStart.height + deltaY;
                        newX = resizeStart.nodeX + deltaX;
                        break;
                    case "bottom-right":
                        newWidth = resizeStart.width + deltaX;
                        newHeight = resizeStart.height + deltaY;
                        break;
                }

                element.style.width = `${Math.max(100, newWidth)}px`;
                element.style.height = `${Math.max(60, newHeight)}px`;
                if (onPositionChange) {
                    onPositionChange(node.uid, newX, newY);
                }
            }
        };

        const handleMouseUp = () => {
            setIsDragging(false);
            setIsResizing(null);
        };

        document.addEventListener("mousemove", handleMouseMove);
        document.addEventListener("mouseup", handleMouseUp);

        return () => {
            document.removeEventListener("mousemove", handleMouseMove);
            document.removeEventListener("mouseup", handleMouseUp);
        };
    }, [isDragging, isResizing, dragStart, resizeStart, node.uid, onPositionChange]);

    const handleMouseDown = useCallback((e: MouseEvent<HTMLDivElement>) => {
        if (e.button !== 0) return;

        setDragStart({ x: e.clientX - node.x, y: e.clientY - node.y });
        setIsDragging(true);
    }, [node.x, node.y]);

    const handleCornerMouseDown = useCallback((corner: ResizeCorner, e: MouseEvent<HTMLDivElement>) => {
        e.stopPropagation();
        const element = nodeRef.current;
        if (!element) return;

        const rect = element.getBoundingClientRect();
        setResizeStart({
            x: e.clientX,
            y: e.clientY,
            width: rect.width,
            height: rect.height,
            nodeX: node.x,
            nodeY: node.y
        });
        setIsResizing(corner);
    }, [node.x, node.y]);

    const cornerSize = 12;
    const cornerStyle = {
        position: "absolute" as const,
        width: `${cornerSize}px`,
        height: `${cornerSize}px`,
        background: "#4A90E2",
        border: "1px solid white",
        borderRadius: "2px",
        cursor: "pointer",
        zIndex: 10
    };

    return (
        <div
            ref={nodeRef}
            onMouseDown={handleMouseDown}
            style={{
                position: "absolute",
                left: `${node.x}px`,
                top: `${node.y}px`,
                minWidth: "200px",
                minHeight: "60px",
                padding: "16px",
                background: "white",
                border: "1px solid #ccc",
                borderRadius: "8px",
                boxSizing: "border-box",
                cursor: isDragging ? "grabbing" : "grab",
                userSelect: "none"
            }}
        >
            <div
                onMouseDown={(e) => handleCornerMouseDown("top-left", e)}
                style={{ ...cornerStyle, top: `-${cornerSize / 2}px`, left: `-${cornerSize / 2}px`, cursor: "nwse-resize" }}
            />
            <div
                onMouseDown={(e) => handleCornerMouseDown("top-right", e)}
                style={{ ...cornerStyle, top: `-${cornerSize / 2}px`, right: `-${cornerSize / 2}px`, cursor: "nesw-resize" }}
            />
            <div
                onMouseDown={(e) => handleCornerMouseDown("bottom-left", e)}
                style={{ ...cornerStyle, bottom: `-${cornerSize / 2}px`, left: `-${cornerSize / 2}px`, cursor: "nesw-resize" }}
            />
            <div
                onMouseDown={(e) => handleCornerMouseDown("bottom-right", e)}
                style={{ ...cornerStyle, bottom: `-${cornerSize / 2}px`, right: `-${cornerSize / 2}px`, cursor: "nwse-resize" }}
            />

            <h3 style={{ margin: "0 0 8px 0", fontSize: "14px" }}>{node.name}</h3>
            <p style={{ margin: 0, fontSize: "12px", color: "#666" }}>{node.description}</p>
        </div>
    );
}