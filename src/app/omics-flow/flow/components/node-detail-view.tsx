import { useEffect, useRef, useState } from "react";
import { Node, RunConfiguration, RunEnvironment } from "../components/node";
import RunnerTerminal, { RunnerTerminalRef } from "./runner-terminal";
import { useSshCredentials } from "../../hook/use-ssh-credentials";
import FileBrowser from "./file-browser";

export default function NodeDetailView({ node }: { node: Node }) {
    const terminalRef = useRef<RunnerTerminalRef>(null);
    const credentials = useSshCredentials();
    const [environments, setEnvironments] = useState<RunEnvironment[]>([]);
    const [environmentIsLoading, setEnvironmentIsLoading] = useState(false);
    const [selectedEnvironment, setSelectedEnvironment] = useState<string>('');
    const [pendingRunConfig, setPendingRunConfig] = useState<RunConfiguration | null>(null);
    const [pendingArgs, setPendingArgs] = useState<string>('');
    const [name, setName] = useState(node.name);
    const [description, setDescription] = useState(node.description);

    useEffect(() => {
        console.log(node);
    });

    useEffect(() => {
        const timer = setTimeout(() => {
            if (name !== node.name) {
                updateNodeField('name', name);
            }
        }, 500);
        return () => clearTimeout(timer);
    }, [name]);

    useEffect(() => {
        const timer = setTimeout(() => {
            if (description !== node.description) {
                updateNodeField('description', description);
            }
        }, 500);
        return () => clearTimeout(timer);
    }, [description]);

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
    }, []);

    const handleRunBash = () => {
        if (!terminalRef.current) return;

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: selectedEnvironment || 'python-scimap/singularity',
        });
    };

    const handleRunFromFileBrowser = (config: RunConfiguration) => {
        setPendingRunConfig(config);
        setPendingArgs(config.args.join('\n'));
    };

    const handleExecuteScript = async () => {
        if (!terminalRef.current || !pendingRunConfig) return;

        const parsedArgs = pendingArgs ? pendingArgs.split('\n').filter(arg => arg.trim()) : [];

        const updatedConfig: RunConfiguration = {
            ...pendingRunConfig,
            args: parsedArgs
        };

        await handleRunConfigChange(pendingRunConfig.scriptPath, updatedConfig);

        terminalRef.current.run({
            node_id: node.uid,
            linux_user: credentials.linux_user,
            private_key: credentials.private_key,
            environment: updatedConfig.environment,
            script: updatedConfig.scriptPath,
            args: updatedConfig.args,
        });

        setPendingRunConfig(null);
        setPendingArgs('');
    };

    const handleCancelRun = () => {
        setPendingRunConfig(null);
        setPendingArgs('');
    };

    const handleRunConfigChange = async (scriptPath: string, config: RunConfiguration) => {
        const updatedConfigs = {
            ...(node.runConfigurations || {}),
            [scriptPath]: config
        };

        await updateNodeField('run-configurations', updatedConfigs);
    };

    return (
        <div style={{ padding: "20px" }}>
            <textarea
                value={name}
                onChange={(e) => setName(e.target.value)}
                style={{
                    width: "100%",
                    fontSize: "1.5em",
                    fontWeight: "bold",
                    border: "1px solid #ccc",
                    borderRadius: "4px",
                    padding: "8px",
                    marginBottom: "10px",
                    fontFamily: "inherit",
                    resize: "vertical",
                    minHeight: "40px"
                }}
            />
            <textarea
                value={description}
                onChange={(e) => setDescription(e.target.value)}
                style={{
                    width: "100%",
                    fontSize: "1em",
                    border: "1px solid #ccc",
                    borderRadius: "4px",
                    padding: "8px",
                    marginBottom: "10px",
                    fontFamily: "inherit",
                    resize: "vertical",
                    minHeight: "60px"
                }}
            />
            <h3>Exports:</h3>
            <FileBrowser
                node_id={node.uid}
                linux_user={credentials.linux_user}
                private_key={credentials.private_key}
                run_configurations={node.runConfigurations || {}}
                environments={environments}
                onRun={handleRunFromFileBrowser}
                onRunConfigChange={handleRunConfigChange}
                environment='python-scimap'
            />

            <hr/>
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

                <button
                    onClick={handleRunBash}
                    disabled={terminalRef.current?.isRunning}
                >
                    {terminalRef.current?.isRunning ? 'Running...' : 'Run Bash'}
                </button>
            </div>
            <hr/>

            {pendingRunConfig && (
                <div style={{ marginTop: "20px", padding: "10px", border: "1px solid #4A90E2" }}>
                    <div><strong>Script:</strong> {pendingRunConfig.scriptPath} | <strong>Env:</strong> {pendingRunConfig.environment}</div>
                    <textarea
                        value={pendingArgs}
                        onChange={(e) => setPendingArgs(e.target.value)}
                        rows={3}
                        style={{ width: "100%", fontFamily: "monospace", marginTop: "5px" }}
                        placeholder="--input-selector ../001-raw-input/output/*/*.csv --output ./output/cell-counts.csv"
                    />
                    <button onClick={handleExecuteScript} disabled={terminalRef.current?.isRunning}>
                        Run Script
                    </button>
                    <button onClick={handleCancelRun} style={{ marginLeft: "5px" }}>
                        Cancel
                    </button>
                </div>
            )}

            <RunnerTerminal ref={terminalRef} />
        </div>
    )

}