import { useRef } from "react";
import { Node } from "../components/node";
import RunnerTerminal, { RunnerTerminalRef } from "./runner-terminal";
import { useSshCredentials } from "../../hook/use-ssh-credentials";

export default function NodeDetailView({ node }: { node: Node }) {
    const terminalRef = useRef<RunnerTerminalRef>(null);
    const credentials = useSshCredentials();

    const handleRunScript = () => {
        if (!terminalRef.current) return;

        // TODO: Replace these with actual properties from the node object
        const scriptPath = 'compute-stats.py'; // Where is this in the node?
        const scriptArgs = [
            "--input-selector", "../001-raw-input/output/*/*.csv",
            "--output", "/nodes/002-cell-stats/output/cell-counts.csv"
        ]; // Where are these in the node?

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: 'python-scimap',
            script: scriptPath,
            args: scriptArgs,
        });
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
            <button
                onClick={handleRunScript}
                disabled={terminalRef.current?.isRunning}
            >
                {terminalRef.current?.isRunning ? 'Running...' : 'Run Script'}
            </button>

            <RunnerTerminal ref={terminalRef} />
        </div>
    )

}