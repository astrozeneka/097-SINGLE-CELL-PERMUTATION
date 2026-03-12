"use client";

import { DrawingCanvas } from "./drawing-canvas";
import { UnderlyingCanvas } from "./underlying-canvas";

export function SampleCutterClient({ data }: { data: number[] }) {
    return (
        <div>
            <UnderlyingCanvas data={data} />
            <DrawingCanvas />
        </div>
    );
}
