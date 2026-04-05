import { useState } from "react";
import { Node } from "../components/node";

export default function NodeDetailView({ node }: { node: Node }) {
    
    const [isScriptRunning, setIsScriptRunning] = useState(false);
    const handleRunScript = () => {

    };

    return (

        <div style={{ padding: "20px" }}>
            <h2>{node.name}</h2>
            <p>{node.description}</p>
            <h3>Exports:</h3>
            <ul>
                {Object.entries(node.exports).map(([key, value]) => (
                    <li key={key}>{key}: {value}</li>
                ))}
            </ul>
            <button onClick={handleRunScript}>Run Script</button>
        </div>
    )

}