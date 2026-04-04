import { Node } from "./node";

interface EdgeProps {
    source: Node;
    destination: Node;
}

export function Edge({ source, destination }: EdgeProps) {
    const nodeWidth = 200;
    const nodeHeight = 60;

    const startX = source.x + nodeWidth;
    const startY = source.y + nodeHeight / 2;
    const endX = destination.x;
    const endY = destination.y + nodeHeight / 2;

    const midX = (startX + endX) / 2;

    const path = `M ${startX} ${startY} C ${midX} ${startY}, ${midX} ${endY}, ${endX} ${endY}`;

    return (
        <svg
            style={{
                position: "absolute",
                top: 0,
                left: 0,
                width: "100%",
                height: "100%",
                pointerEvents: "none",
                overflow: "visible"
            }}
        >
            <defs>
                <marker
                    id={`arrowhead-${source.uid}-${destination.uid}`}
                    markerWidth="10"
                    markerHeight="10"
                    refX="9"
                    refY="3"
                    orient="auto"
                >
                    <polygon points="0 0, 10 3, 0 6" fill="#999" />
                </marker>
            </defs>
            <path
                d={path}
                stroke="#999"
                strokeWidth="2"
                fill="none"
                markerEnd={`url(#arrowhead-${source.uid}-${destination.uid})`}
            />
        </svg>
    );
}
