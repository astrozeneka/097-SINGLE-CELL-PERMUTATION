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
    const [selectedEnvironment, setSelectedEnvironment] = useState<string>('');
    const [scriptArgs, setScriptArgs] = useState<string>('');
    const debounceTimerRef = useRef<NodeJS.Timeout | null>(null);

    const updateNodeField = async (field: string, value: any) => {
        try {
            await fetch(`http://192.168.64.3:3000/nodes/${node.uid}`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    linux_user: credentials.linux_user,
                    private_key: credentials.private_key,
                    [field]: value
                })
            });
        } catch (error) {
            console.error(`Error updating node: ${error}`);
        }
    };

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
        setScriptArgs((node as any).scriptArgs || '');
    }, []);

    const handleScriptArgsChange = (value: string) => {
        setScriptArgs(value);

        if (debounceTimerRef.current) {
            clearTimeout(debounceTimerRef.current);
        }

        debounceTimerRef.current = setTimeout(() => {
            updateNodeField('scriptArgs', value);
        }, 500);
    };

    const handleRunBash = () => {
        if (!terminalRef.current) return;

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: selectedEnvironment || 'python-scimap/singularity',
        });
    };

    const handleRunScript = () => {
        if (!terminalRef.current) return;

        const scriptPath = 'compute-stats.py';
        const parsedArgs = scriptArgs ? scriptArgs.split('\n').filter(arg => arg.trim()) : [];

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: selectedEnvironment || 'python-scimap/singularity',
            script: scriptPath,
            args: parsedArgs,
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
            <div style={{ marginBottom: "10px" }}>
                <label htmlFor="environment-selector">Environment: </label>
                <select
                    id="environment-selector"
                    value={selectedEnvironment}
                    onChange={(e) => setSelectedEnvironment(e.target.value)}
                >
                    <option value="">Select environment</option>
                    {environments.map(env => (
                        <option key={env.name} value={env.name}>
                            {env.name} ({env.type})
                        </option>
                    ))}
                </select>
            </div>
            <div style={{ marginBottom: "10px" }}>
                <label htmlFor="script-args">Script Arguments (one per line):</label>
                <br />
                <textarea
                    id="script-args"
                    value={scriptArgs}
                    onKeyUp={(e) => handleScriptArgsChange(e.currentTarget.value)}
                    onChange={(e) => setScriptArgs(e.target.value)}
                    rows={5}
                    style={{ width: "100%", fontFamily: "monospace" }}
                    placeholder="--input-selector&#10;../001-raw-input/output/*/*.csv&#10;--output&#10;./output/cell-counts.csv"
                />
            </div>
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