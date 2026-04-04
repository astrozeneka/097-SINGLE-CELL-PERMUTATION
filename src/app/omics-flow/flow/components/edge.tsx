import { Node, NodeDimensions } from "./node";

interface EdgeProps {
    source: Node;
    destination: Node;
    sourceDimensions?: NodeDimensions;
    destinationDimensions?: NodeDimensions;
    startOffset?: number;
    endOffset?: number;
}

type Side = "top" | "right" | "bottom" | "left";

interface ConnectionPoint {
    x: number;
    y: number;
    side: Side;
}

interface RouteResult {
    start: ConnectionPoint;
    end: ConnectionPoint;
    distance: number;
}

interface PathData {
    path: string;
    svgBounds: {
        x: number;
        y: number;
        width: number;
        height: number;
    };
}

function getConnectionPoints(node: Node, dimensions: NodeDimensions): Record<Side, ConnectionPoint> {
    const cx = node.x + dimensions.width / 2;
    const cy = node.y + dimensions.height / 2;

    return {
        top: { x: cx, y: node.y, side: "top" },
        right: { x: node.x + dimensions.width, y: cy, side: "right" },
        bottom: { x: cx, y: node.y + dimensions.height, side: "bottom" },
        left: { x: node.x, y: cy, side: "left" }
    };
}

function calculateDistance(p1: ConnectionPoint, p2: ConnectionPoint): number {
    return Math.sqrt(Math.pow(p2.x - p1.x, 2) + Math.pow(p2.y - p1.y, 2));
}

function findBestRoute(
    sourcePoints: Record<Side, ConnectionPoint>,
    destPoints: Record<Side, ConnectionPoint>
): RouteResult {
    let bestRoute: RouteResult | null = null;
    const sides: Side[] = ["top", "right", "bottom", "left"];

    for (const sourceSide of sides) {
        for (const destSide of sides) {
            const start = sourcePoints[sourceSide];
            const end = destPoints[destSide];
            const distance = calculateDistance(start, end);

            if (!bestRoute || distance < bestRoute.distance) {
                bestRoute = { start, end, distance };
            }
        }
    }

    return bestRoute!;
}

function createPathData(
    source: Node,
    destination: Node,
    sourceDimensions: NodeDimensions,
    destinationDimensions: NodeDimensions,
    route: RouteResult,
    startOffset: number,
    endOffset: number
): PathData {
    const { start, end } = route;

    const srcCenterX = source.x + sourceDimensions.width / 2;
    const srcCenterY = source.y + sourceDimensions.height / 2;
    const dstCenterX = destination.x + destinationDimensions.width / 2;
    const dstCenterY = destination.y + destinationDimensions.height / 2;

    const svgX = Math.min(srcCenterX, dstCenterX);
    const svgY = Math.min(srcCenterY, dstCenterY);
    const svgWidth = Math.abs(dstCenterX - srcCenterX);
    const svgHeight = Math.abs(dstCenterY - srcCenterY);

    const toLocal = (x: number, y: number) => ({
        x: x - svgX,
        y: y - svgY
    });

    let exitX = start.x;
    let exitY = start.y;

    switch (start.side) {
        case "right":
            exitX += startOffset;
            break;
        case "left":
            exitX -= startOffset;
            break;
        case "bottom":
            exitY += startOffset;
            break;
        case "top":
            exitY -= startOffset;
            break;
    }

    let entryX = end.x;
    let entryY = end.y;

    switch (end.side) {
        case "right":
            entryX += endOffset;
            break;
        case "left":
            entryX -= endOffset;
            break;
        case "bottom":
            entryY += endOffset;
            break;
        case "top":
            entryY -= endOffset;
            break;
    }

    const controlOffset = Math.min(Math.abs(entryX - exitX), Math.abs(entryY - exitY)) / 2;

    let cp1X = exitX;
    let cp1Y = exitY;
    let cp2X = entryX;
    let cp2Y = entryY;

    if (start.side === "right" || start.side === "left") {
        cp1X += (start.side === "right" ? 1 : -1) * controlOffset;
    } else {
        cp1Y += (start.side === "bottom" ? 1 : -1) * controlOffset;
    }

    if (end.side === "right" || end.side === "left") {
        cp2X += (end.side === "right" ? 1 : -1) * controlOffset;
    } else {
        cp2Y += (end.side === "bottom" ? 1 : -1) * controlOffset;
    }

    const srcCenter = toLocal(srcCenterX, srcCenterY);
    const exit = toLocal(exitX, exitY);
    const cp1 = toLocal(cp1X, cp1Y);
    const cp2 = toLocal(cp2X, cp2Y);
    const entry = toLocal(entryX, entryY);
    const dstCenter = toLocal(dstCenterX, dstCenterY);

    const path = `M ${srcCenter.x} ${srcCenter.y} L ${exit.x} ${exit.y} C ${cp1.x} ${cp1.y}, ${cp2.x} ${cp2.y}, ${entry.x} ${entry.y} L ${dstCenter.x} ${dstCenter.y}`;

    return {
        path,
        svgBounds: { x: svgX, y: svgY, width: svgWidth, height: svgHeight }
    };
}

export function Edge({
    source,
    destination,
    sourceDimensions,
    destinationDimensions,
    startOffset = 10,
    endOffset = 20
}: EdgeProps) {
    if (!sourceDimensions || !destinationDimensions) {
        return null;
    }

    const sourcePoints = getConnectionPoints(source, sourceDimensions);
    const destPoints = getConnectionPoints(destination, destinationDimensions);
    const route = findBestRoute(sourcePoints, destPoints);
    const { path, svgBounds } = createPathData(
        source,
        destination,
        sourceDimensions,
        destinationDimensions,
        route,
        startOffset,
        endOffset
    );

    return (
        <svg
            style={{
                position: "absolute",
                left: `${svgBounds.x}px`,
                top: `${svgBounds.y}px`,
                width: `${svgBounds.width}px`,
                height: `${svgBounds.height}px`,
                pointerEvents: "none",
                overflow: "visible",
                border: "1px solid pink"
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
