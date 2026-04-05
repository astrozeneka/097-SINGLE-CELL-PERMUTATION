"use client";

import { useState } from "react";
import { Grid } from "./components/grid";
import { Node } from "./components/node";
import NodeDetailView from "./components/node-detail-view";

export default function Flow() {
    // The node that is currently selected by the user
    const [selectedNode, setSelectedNode] = useState<Node | null>(null);

    return (
        <div style={{display: "flex", flexDirection: "row", width: "100vw", height: "100vh", position: "absolute", top: 0, left: 0}}>
            <div style={{ flexGrow: 1, height: "100%"}}>
                <Grid selectedNode={selectedNode} setSelectedNode={setSelectedNode}></Grid>
            </div>
            <div style={{ flex: 1}}>
                {selectedNode ? (
                    <NodeDetailView node={selectedNode} />
                ) : (
                    <div style={{ padding: "20px" }}>
                        <h2>No node selected</h2>
                        <p>Please select a node to see its details.</p>
                    </div>
                )}
            </div>
        </div>
    )
}