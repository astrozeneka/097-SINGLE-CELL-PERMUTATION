export interface NodeInfo {
    uid: string;
    name: string;
    description: string;
    "canvas-x": number;
    "canvas-y": number;
    "edges-in": string[];
    exports: Record<string, string>;
}

export interface NodeData {
    info: NodeInfo;
}

export class Node {
    uid: string;
    name: string;
    description: string;
    x: number;
    y: number;
    edgesIn: string[];
    exports: Record<string, string>;

    constructor(data: NodeData) {
        this.uid = data.info.uid;
        this.name = data.info.name;
        this.description = data.info.description;
        this.x = data.info["canvas-x"];
        this.y = data.info["canvas-y"];
        this.edgesIn = data.info["edges-in"];
        this.exports = data.info.exports;
    }
}

export function NodeComponent({ node }: { node: Node }) {
    return (
        <div
            style={{
                position: "absolute",
                left: `${node.x}px`,
                top: `${node.y}px`,
                padding: "16px",
                background: "white",
                border: "1px solid #ccc",
                borderRadius: "8px",
                minWidth: "200px",
            }}
        >
            <h3 style={{ margin: "0 0 8px 0" }}>{node.name}</h3>
            <p style={{ margin: 0, fontSize: "14px", color: "#666" }}>{node.description}</p>
        </div>
    );
}