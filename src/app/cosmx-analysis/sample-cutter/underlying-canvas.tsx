"use client";

import React from "react";

// HARD-CODING: FOR TESTING ONLY
import LarcBSampleDataCSV from "./larc-b-sample-data-csv";
const arr = LarcBSampleDataCSV.split("\n").slice(1).map(line => {
    const [id, x, y, cluster] = line.split(",");
    return { id, x: parseFloat(x), y: parseFloat(y), cluster: parseInt(cluster) };
});

export function UnderlyingCanvas() {
    const canvasRef = React.useRef<HTMLCanvasElement>(null);
    return (
        <canvas ref={canvasRef} width={500} height={500} />
    )
}