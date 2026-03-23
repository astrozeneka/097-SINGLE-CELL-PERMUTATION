
interface PolygonManagerCanvasProps {
    size: { w: number; h: number };
    transform: any;
    onTransform: (transform: any) => void;
}

const BRUSH_RADIUS = 100;
const BRUSH_VERTS = 16;

export default function PolygonManagerCanvas({ size, transform, onTransform }: PolygonManagerCanvasProps) {
    return (
        <canvas style={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%" }} />
    )
}