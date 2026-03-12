import React from "react";

export function DrawingCanvas() {
    const canvasRef = React.useRef<HTMLCanvasElement>(null);
    return (
        <canvas ref={canvasRef} width={500} height={500} />
    )
}