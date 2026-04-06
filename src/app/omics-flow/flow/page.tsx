"use client";

import { useState } from "react";
import { Grid } from "./components/grid";
import { Node } from "./components/node";
import NodeDetailView from "./components/node-detail-view";
import HorizontalSplit from "./components/horizontal-split";

export default function Flow() {
    // The nodes that are currently selected by the user
    const [selectedNodes, setSelectedNodes] = useState<Node[]>([]);

    return (
        <div style={{ display: "inline-block", position: "absolute", width: "100vw", height: "100vh" }}>
            <HorizontalSplit>
                <div style={{ flexGrow: 1, height: "100%"}}>
                    <Grid selectedNodes={selectedNodes} setSelectedNodes={setSelectedNodes}></Grid>
                </div>
                <div style={{ flex: 1}}>
                    {selectedNodes.length === 1 ? (
                        <NodeDetailView key={selectedNodes[0].uid} node={selectedNodes[0]} />
                    ) : (selectedNodes.length === 0) ? (
                        <div style={{ padding: "20px" }}>
                            <h2>No node selected</h2>
                            <p>Please select a node to see its details.</p>
                        </div>
                    ) : (
                        <div style={{ padding: "20px" }}>
                            <h2>{selectedNodes.length} nodes selected</h2>
                        </div>
                    )}
                </div>
            </HorizontalSplit>
        </div>
    )
}