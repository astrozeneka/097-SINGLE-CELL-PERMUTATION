"use client";

import { useState } from "react";
import { Grid } from "./components/grid";
import { Node } from "./components/node";
import NodeDetailView from "./components/node-detail-view";
import HorizontalSplit from "./components/horizontal-split";
import { CreateNodeDialog } from "./components/create-node-dialog";
import { NodeManager } from "./providers/node-manager";
import { useSshCredentials } from "../hook/use-ssh-credentials";

export default function Flow() {
    // The nodes that are currently selected by the user
    const [selectedNodes, setSelectedNodes] = useState<Node[]>([]);
    const [isFloatingDialogOpen, setIsFloatingDialogOpen] = useState(false);
    const [, forceUpdate] = useState({});
    const credentials = useSshCredentials();
    const nodeManager = NodeManager.instance;

    const handleFloatingCreate = async (name: string, description: string) => {
        await nodeManager.createNode(
            name,
            description,
            0,
            0,
            credentials.linux_user,
            credentials.private_key
        );
        forceUpdate({});
    };

    return (
        <div style={{ display: "inline-block", position: "absolute", width: "100vw", height: "100vh" }}>
            <HorizontalSplit>
                <div style={{ flexGrow: 1, height: "100%"}}>
                    <Grid selectedNodes={selectedNodes} setSelectedNodes={setSelectedNodes}></Grid>
                </div>
                <div style={{ flex: 1, height: "100vh", overflow: "scroll"}}>
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

            {/* Floating button */}
            <button
                onClick={() => setIsFloatingDialogOpen(true)}
                style={{
                    position: "fixed",
                    bottom: "24px",
                    right: "24px",
                    width: "56px",
                    height: "56px",
                    borderRadius: "50%",
                    background: "#4A90E2",
                    border: "none",
                    color: "white",
                    fontSize: "24px",
                    cursor: "pointer",
                    boxShadow: "0 4px 12px rgba(0, 0, 0, 0.15)",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    zIndex: 1000
                }}
                title="Create new node"
            >
                +
            </button>

            <CreateNodeDialog
                isOpen={isFloatingDialogOpen}
                onClose={() => setIsFloatingDialogOpen(false)}
                onCreate={handleFloatingCreate}
            />
        </div>
    )
}