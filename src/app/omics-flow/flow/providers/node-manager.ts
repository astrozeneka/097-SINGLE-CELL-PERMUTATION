import { Node, NodeData } from "../components/node";

const HARDCODED_NODE_DATA: NodeData[] = [
    {
        "dir": "nodes/001-raw-input",
        "info": {
            "uid": "001-raw-input",
            "name": "Raw input - Sarcoma",
            "description": "Raw input data",
            "canvas-x": 100,
            "canvas-y": 100,
            "edges-in": [],
            "exports": {
                "raw": "*/*.csv"
            }
        }
    },
    {
        "dir": "nodes/002-processing",
        "info": {
            "uid": "002-processing",
            "name": "Data Processing",
            "description": "Process raw data",
            "canvas-x": 400,
            "canvas-y": 150,
            "edges-in": ["001-raw-input"],
            "exports": {
                "processed": "*/*.json"
            }
        }
    },
    {
        "dir": "nodes/003-analysis",
        "info": {
            "uid": "003-analysis",
            "name": "Data Analysis",
            "description": "Analyze processed data",
            "canvas-x": 700,
            "canvas-y": 200,
            "edges-in": ["002-processing"],
            "exports": {
                "results": "*/*.txt"
            }
        }
    }
];

export class NodeManager {
    private nodes: Node[] = [];

    private constructor() {
        this.nodes = HARDCODED_NODE_DATA.map(data => new Node(data));
    }

    public static readonly instance: NodeManager = new NodeManager();

    public getNodes(): Node[] {
        return this.nodes;
    }
}