import { useEffect, useRef, useState } from "react";
import { Node } from "../components/node";
import RunnerTerminal, { RunnerTerminalRef } from "./runner-terminal";
import { useSshCredentials } from "../../hook/use-ssh-credentials";
import FileBrowser from "./file-browser";

interface RunEnvironment {
    name: string,
    type: "docker" | "singularity" | "conda" | "system",
    status: "up-to-date" | "needs-rebuild" | "no-build"
}

export default function NodeDetailView({ node }: { node: Node }) {
    const terminalRef = useRef<RunnerTerminalRef>(null);
    const credentials = useSshCredentials();
    const [environments, setEnvironments] = useState<RunEnvironment[]>([]);
    const [environmentIsLoading, setEnvironmentIsLoading] = useState(false);

    useEffect(() => {
        const fetchEnvironments = async () => {
            console.log("Fetch environments")
            // Load environment list from the backend
            setEnvironmentIsLoading(true);
            try {
                const params = new URLSearchParams({
                    linux_user: credentials.linux_user,
                    private_key: credentials.private_key,
                });
                const response = await fetch(`http://192.168.64.3:3000/environments?${params}`, {
                    method: 'GET',
                    headers: {
                        'Content-Type': 'application/json'
                    }
                })
                if (!response.ok) {
                    console.error(`Failed to fetch environments: ${response.statusText}`);
                    return;
                }
                const data: { environments: RunEnvironment[] } = await response.json();
                setEnvironments(data.environments);
            } catch (error) {
                console.error(`Error fetching environments: ${error}`);
            } finally {
                setEnvironmentIsLoading(false);
            }
        };

        fetchEnvironments();
    }, []);

    const handleRunBash = () => {
        if (!terminalRef.current) return;

        // TODO: Replace these with actual properties from the node object

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: 'python-scimap/singularity',
        });
    };

    const handleRunScript = () => {
        if (!terminalRef.current) return;


        const scriptPath = 'compute-stats.py'; // Where is this in the node?

        const scriptArgs = [
            "--input-selector", "../001-raw-input/output/*/*.csv",
            "--output", "./output/cell-counts.csv"
        ]; // Where are these in the node?
        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: 'python-scimap/singularity',
            script: scriptPath,
            args: scriptArgs,
        });
    }

    return (
        <div style={{ padding: "20px" }}>
            <h2>{node.name}</h2>
            <p>{node.description}</p>
            <h3>Exports:</h3>
            <FileBrowser
                node_id={node.uid}
                linux_user={credentials.linux_user}
                private_key={credentials.private_key}
                environment='python-scimap'
            />
            <div>
                {environments.length === 0 ? (
                    <p>No environments found.</p>
                ) : (
                    <ul>
                        {environments.map(env => (
                            <li key={env.name}>
                                {env.name} ({env.type}) - {env.status}
                            </li>
                        ))}
                    </ul>
                )}
            </div>
            <ul>
                {Object.entries(node.exports).map(([key, value]) => (
                    <li key={key}>{key}: {value}</li>
                ))}
            </ul>
            <button
                onClick={handleRunBash}
                disabled={terminalRef.current?.isRunning}
            >
                {terminalRef.current?.isRunning ? 'Running...' : 'Run Bash'}
            </button>
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