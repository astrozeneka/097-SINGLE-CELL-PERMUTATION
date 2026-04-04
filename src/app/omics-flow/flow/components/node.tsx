
export class Node {
    title: string;
    description: string;

    constructor(title: string, description: string) {
        this.title = title;
        this.description = description;
    }
}

function NodeComponent({ node }: { node: Node }) {
    return (
        <div>
            <h3>{node.title}</h3>
            <p>{node.description}</p>
        </div>
    );
}