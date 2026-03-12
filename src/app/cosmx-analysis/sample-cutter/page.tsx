import React from "react";

class Point {
    constructor(public x: number, public y: number) {}
}

class Polygon {
    constructor(public vertices: Point[]) {}   
}

export default function SampleCutterPage() {
    const canvasRef = React.useRef<HTMLCanvasElement>(null);
    const currentPolygon = React.useRef<Polygon | null>(null);

    React.useEffect(() => {
        const canvas = canvasRef.current;
        if (canvas) {
            const context = canvas.getContext('2d');
            if (context) {
                context.fillStyle = 'red';
                context.fillRect(0, 0, canvas.width, canvas.height);
            }
        }
    }, []);

    const handleClick = () => {
        // Begin a polygon if not already started,
        // Or add points
    };

    return (
        <div>
            <canvas ref={canvasRef} onClick={handleClick} width="400" height="400" style={{ border: '1px solid black' }}></canvas>
            <div>
                {/* A ready to copy geojson will be here */}
            </div>
        </div>
    );
}