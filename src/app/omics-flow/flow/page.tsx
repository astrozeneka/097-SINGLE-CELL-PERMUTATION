"use client";

import { useState } from "react";
import { Grid } from "./components/grid";
import { Node } from "./components/node";

export default function Flow() {
    // The node that is currently selected by the user
    const [selectedNode, setSelectedNode] = useState<Node | null>(null);

    return (
        <div style={{display: "flex", flexDirection: "row", width: "100vw", height: "100vh", position: "absolute", top: 0, left: 0}}>
            <div style={{ flexGrow: 1, height: "100%"}}>
                <Grid selectedNode={selectedNode} setSelectedNode={setSelectedNode}></Grid>
            </div>
            <div style={{ flex: 1}}>
                Hello
            </div>
        </div>
    )
}