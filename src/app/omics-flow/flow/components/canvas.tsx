import { useRef, useState, useCallback, useEffect, MouseEvent, ReactNode } from "react";

export interface Transform {
    x: number;
    y: number;
    scale: number;
}

interface CanvasProps {
    children: ReactNode;
    transform: Transform;
    setTransform: (transform: Transform | ((prev: Transform) => Transform)) => void;
}

export function Canvas({ children, transform, setTransform }: CanvasProps) {
    const canvasRef = useRef<HTMLDivElement>(null);
    const [isPanning, setIsPanning] = useState(false);
    const [panStart, setPanStart] = useState({ x: 0, y: 0 });

    useEffect(() => {
        const canvas = canvasRef.current;
        if (!canvas) return;

        const handleWheel = (e: globalThis.WheelEvent) => {
            e.preventDefault();
            e.stopPropagation();

            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;

            // Calculate zoom
            const delta = e.deltaY;
            const zoomIntensity = 0.001;
            const scaleFactor = 1 - delta * zoomIntensity;

            setTransform(prev => {
                const newScale = Math.min(Math.max(0.1, prev.scale * scaleFactor), 5);

                // Zoom towards mouse position
                const scaleRatio = newScale / prev.scale;
                const newX = mouseX - (mouseX - prev.x) * scaleRatio;
                const newY = mouseY - (mouseY - prev.y) * scaleRatio;

                return {
                    x: newX,
                    y: newY,
                    scale: newScale
                };
            });
        };

        canvas.addEventListener('wheel', handleWheel, { passive: false });

        return () => {
            canvas.removeEventListener('wheel', handleWheel);
        };
    }, [setTransform]);

    const handleMouseDown = useCallback((e: MouseEvent<HTMLDivElement>) => {
        // Right mouse button (button 2) for panning
        if (e.button === 2) {
            e.preventDefault();
            setIsPanning(true);
            setPanStart({ x: e.clientX - transform.x, y: e.clientY - transform.y });
        }
    }, [transform.x, transform.y]);

    const handleMouseMove = useCallback((e: MouseEvent<HTMLDivElement>) => {
        if (!isPanning) return;

        setTransform(prev => ({
            ...prev,
            x: e.clientX - panStart.x,
            y: e.clientY - panStart.y
        }));
    }, [isPanning, panStart]);

    const handleMouseUp = useCallback((e: MouseEvent<HTMLDivElement>) => {
        if (e.button === 2) {
            setIsPanning(false);
        }
    }, []);

    return (
        <div
            ref={canvasRef}
            onMouseDown={handleMouseDown}
            onMouseMove={handleMouseMove}
            onMouseUp={handleMouseUp}
            onMouseLeave={() => setIsPanning(false)}
            onContextMenu={(e) => e.preventDefault()}
            style={{
                width: "100%",
                height: "100%",
                position: "relative",
                overflow: "hidden",
                cursor: isPanning ? "grabbing" : "default",
                touchAction: "none"
            }}
        >
            <div
                style={{
                    transform: `translate(${transform.x}px, ${transform.y}px) scale(${transform.scale})`,
                    transformOrigin: "0 0",
                    width: "100%",
                    height: "100%",
                    position: "relative",
                }}
            >
                {children}
            </div>
        </div>
    );
}
