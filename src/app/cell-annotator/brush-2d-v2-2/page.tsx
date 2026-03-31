"use client";

import HorizontalSplit from "../horizontal-split";

export default function Page() {
    return (
        <div style={{ width: "100vw", height: "100vh", position: "absolute", left: 0, top: 0 }}>
            <HorizontalSplit>
                <div style={{ background: "lightblue" }}>Hello world left</div>
                <div style={{ background: "lightcoral" }}>Hello world right</div>
            </HorizontalSplit>
        </div>
    );
}