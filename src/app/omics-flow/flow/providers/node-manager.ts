import { Node, NodeData } from "../components/node";
import { OMICS_FLOW_API_URL } from "@/config/api";


export class NodeManager {
    // private nodes: Node[] = [];
    private nodes: Record<string, Node>;

    private constructor() {
        this.nodes = {};
        /*HARDCODED_NODE_DATA.forEach(data => {
            const node = new Node(data);
            this.nodes[node.uid] = node;
        });*/
    }

    public static readonly instance: NodeManager = new NodeManager();

    public getNodes(): Node[] {
        return Object.values(this.nodes);
    }

    public setNodes(nodeDataList: NodeData[]): void {
        this.nodes = {};
        nodeDataList.forEach(data => {
            const node = new Node(data);
            this.nodes[node.uid] = node;
        });
    }

    public addNode(nodeData: NodeData): void {
        const node = new Node(nodeData);
        this.nodes[node.uid] = node;
    }

    public generateCloneName(originalUid: string): string {
        let suffix = 1;
        let newUid = `${originalUid}-clone`;

        while (this.nodes[newUid]) {
            suffix++;
            newUid = `${originalUid}-clone-${suffix}`;
        }

        return newUid;
    }

    public async createNode(
        name: string,
        description: string,
        x: number,
        y: number,
        linuxUser: string,
        privateKey: string
    ): Promise<Node | null> {
        console.log(`Creating node "${name}" at (${x}, ${y}) with user ${linuxUser}`);

        const requestBody = {
            linux_user: linuxUser,
            private_key: privateKey,
            name,
            description,
            x,
            y
        };

        const response = await fetch(`${OMICS_FLOW_API_URL}/nodes`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify(requestBody)
        });

        if (!response.ok) {
            throw new Error(`Failed to create node: ${response.statusText}`);
        }

        const data = await response.json();

        const nodeData: NodeData = {
            info: {
                ...data.node.info,
                exports: data.node.info.exports || {}
            }
        };

        this.addNode(nodeData);
        return this.nodes[data.node.info.uid];
    }

    public async cloneNode(
        sourceUid: string,
        x: number,
        y: number,
        linuxUser: string,
        privateKey: string
    ): Promise<Node | null> {
        console.log(`Cloning node ${sourceUid} at (${x}, ${y}) with user ${linuxUser}`);
        const destinationUid = this.generateCloneName(sourceUid);

        const requestBody = {
            linux_user: linuxUser,
            x,
            y,
            src_node: sourceUid,
            dst_node: destinationUid,
            private_key: privateKey
        };

        const response = await fetch(`${OMICS_FLOW_API_URL}/nodes/clone`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify(requestBody)
        });

        if (!response.ok) {
            throw new Error(`Failed to clone node: ${response.statusText}`);
        }

        const data = await response.json();

        if (data.exitCode !== 0) {
            throw new Error(`Clone failed: ${data.stderr || data.stdout}`);
        }

        const nodeData: NodeData = {
            info: {
                ...data.node.info,
                exports: data.node.info.exports || {}
            }
        };

        this.addNode(nodeData);
        return this.nodes[destinationUid];
    }
}