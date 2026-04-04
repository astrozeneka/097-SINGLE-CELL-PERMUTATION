import { Node } from "../components/node";

export class NodeManager {
    private nodes: Node[] = [];
    private constructor() {
        // Hard code to test it
        this.nodes.push(new Node("Node 1", "This is the first node"));
        this.nodes.push(new Node("Node 2", "This is the second node"));
    }

    public static readonly instance: NodeManager = new NodeManager();
}