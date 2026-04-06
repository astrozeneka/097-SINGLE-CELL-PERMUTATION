'use client';

import { useEffect, useState } from 'react';
import { useFileBrowser } from '@/hooks/use-file-browser';
import { FileItem } from '../lib/file-browser-session';
import { RunConfiguration, RunEnvironment } from './node';

export interface FileBrowserProps {
    node_id: string;
    linux_user: string;
    private_key: string;
    host?: string;
    port?: number;
    wsUrl?: string;
    initialPath?: string;
    environment?: string;
    run_configurations: Record<string, RunConfiguration>;
    environments: RunEnvironment[];
    onRun: (config: RunConfiguration) => void;
    onRunConfigChange: (scriptPath: string, config: RunConfiguration) => void;
}

export default function FileBrowser({
    node_id,
    linux_user,
    private_key,
    host,
    port,
    wsUrl,
    initialPath,
    environment,
    run_configurations,
    environments,
    onRun,
    onRunConfigChange
}: FileBrowserProps) {
    const { connect, toggleDirectory, isConnected, fileItems, error } = useFileBrowser();
    const [scriptEnvironments, setScriptEnvironments] = useState<Record<string, string>>({});

    useEffect(() => {
        const envMap: Record<string, string> = {};
        Object.values(run_configurations).forEach(config => {
            envMap[config.scriptPath] = config.environment;
        });
        setScriptEnvironments(envMap);
    }, [run_configurations]);

    useEffect(() => {
        connect({
            node_id,
            linux_user,
            private_key,
            host,
            port,
            wsUrl,
            initialPath,
            environment
        });
    }, [node_id, linux_user, private_key, host, port, wsUrl, initialPath, environment, connect]);

    const handleToggle = (item: FileItem) => {
        if (item.type === 'directory') {
            toggleDirectory(item);
        }
    };

    const handleDownload = (item: FileItem) => {
        console.log('Download placeholder for:', item.fullPath);
    };

    const isExecutable = (fileName: string) => {
        return fileName.endsWith('.sh') || fileName.endsWith('.py') || fileName.endsWith('.R');
    };

    const handleEnvironmentChange = (scriptPath: string, selectedEnv: string) => {
        setScriptEnvironments(prev => ({
            ...prev,
            [scriptPath]: selectedEnv
        }));

        const config: RunConfiguration = {
            scriptPath,
            environment: selectedEnv,
            args: run_configurations[scriptPath]?.args || []
        };
        onRunConfigChange(scriptPath, config);
    };

    const handleRun = (item: FileItem) => {
        const selectedEnv = scriptEnvironments[item.fullPath];
        if (!selectedEnv) {
            console.error('No environment selected for script:', item.fullPath);
            return;
        }

        const config: RunConfiguration = {
            scriptPath: item.fullPath,
            environment: selectedEnv,
            args: run_configurations[item.fullPath]?.args || []
        };
        onRun(config);
    };

    return (
        <div className="w-full h-full bg-gray-900 text-gray-100 p-4 overflow-auto font-mono text-sm">
            <div className="mb-2 text-xs text-gray-400">
                {isConnected ? '● Connected' : '○ Disconnected'}
                {error && <span className="ml-4 text-red-400">Error: {error}</span>}
            </div>

            <div className="space-y-1">
                {fileItems.map((item, index) => (
                    <div
                        key={`${item.fullPath}-${index}`}
                        className="flex items-center hover:bg-gray-800 px-2 py-1 rounded"
                        style={{ paddingLeft: `${item.depth * 20 + 8}px` }}
                    >
                        {item.type === 'directory' ? (
                            <>
                                <button
                                    onClick={() => handleToggle(item)}
                                    className="mr-2 text-gray-400 hover:text-white w-4 flex-shrink-0"
                                >
                                    {item.isExpanded ? '[-]' : '[+]'}
                                </button>
                                <span className="text-blue-400">📁 {item.name}</span>
                            </>
                        ) : (
                            <>
                                <span className="mr-2 w-4 flex-shrink-0"></span>
                                <span className="text-gray-300 flex-1">📄 {item.name}</span>
                                {isExecutable(item.name) && (
                                    <>
                                        <select
                                            value={scriptEnvironments[item.fullPath] || ''}
                                            onChange={(e) => handleEnvironmentChange(item.fullPath, e.target.value)}
                                            className="ml-2 px-2 py-0.5 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-300"
                                        >
                                            <option value="">Select environment</option>
                                            {environments.map(env => (
                                                <option key={env.name} value={env.name}>
                                                    {env.name} ({env.type})
                                                </option>
                                            ))}
                                        </select>
                                        <button
                                            onClick={() => handleRun(item)}
                                            disabled={!scriptEnvironments[item.fullPath]}
                                            className="ml-2 px-2 py-0.5 text-xs bg-blue-700 hover:bg-blue-600 disabled:bg-gray-600 disabled:cursor-not-allowed rounded text-gray-300"
                                        >
                                            Run
                                        </button>
                                    </>
                                )}
                                <button
                                    onClick={() => handleDownload(item)}
                                    className="ml-4 px-2 py-0.5 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-300"
                                >
                                    Download
                                </button>
                            </>
                        )}
                    </div>
                ))}
            </div>

            {fileItems.length === 0 && isConnected && (
                <div className="text-gray-500 text-center mt-8">Loading...</div>
            )}
        </div>
    );
}