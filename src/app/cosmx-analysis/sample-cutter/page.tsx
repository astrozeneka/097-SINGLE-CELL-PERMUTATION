"use client";

import { DrawingCanvas } from "./drawing-canvas";
import { UnderlyingCanvas } from "./underlying-canvas";






export default function SampleCutterPage() {

    return (
        <div>
            <UnderlyingCanvas />
            <DrawingCanvas />
        </div>
    )
}