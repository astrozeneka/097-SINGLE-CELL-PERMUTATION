import { useRef } from "react";
import Scatter3D from "./scatter-3d";

type GlState = { // check it again claude
    gl: WebGLRenderingContext;
    count: number;
    pixelScaleLoc: WebGLUniformLocation;
    pixelOffsetLoc: WebGLUniformLocation;
    sizeLoc: WebGLUniformLocation;
    polygonBuffer: WebGLBuffer;
};

interface ScatterCanvas3DProps {
    data: T[];
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;
    zAccessor: (d: T) => number;

}

export default function Viewer3D<T>() {
    const canvasRef      = useRef<HTMLCanvasElement>(null);
    const glRef          = useRef<GlState | null>(null);

    return <div style={{ display: "flex", flexDirection: "row", width: "100%", height: "100vh", position: "relative" }}>
        <Scatter3D></Scatter3D>
    </div>;
}