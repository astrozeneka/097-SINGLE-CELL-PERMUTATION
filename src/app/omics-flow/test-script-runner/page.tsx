'use client';

import { useEffect, useRef, useState } from 'react';
import { Terminal } from '@xterm/xterm';
import { FitAddon } from '@xterm/addon-fit';
import '@xterm/xterm/css/xterm.css';
import { useScriptRunner } from '@/hooks/use-script-runner';

export default function TestScriptRunner() {
    const terminalRef = useRef<HTMLDivElement>(null);
    const terminal = useRef<Terminal | null>(null);
    const [hasStarted, setHasStarted] = useState(false);
    const isConnectedRef = useRef(false);
    const isRunningRef = useRef(false);

    const { run, sendInput, isConnected, isRunning } = useScriptRunner({
        onStdout: (data) => {
            terminal.current?.write(data);
        },
        onStderr: (data) => {
            terminal.current?.write(`\x1b[31m${data}\x1b[0m`); // Red color for stderr
        },
        onStatus: (message) => {
            terminal.current?.writeln(`\r\n\x1b[36m[${message}]\x1b[0m\r\n`); // Cyan color for status
        },
        onError: (message) => {
            terminal.current?.writeln(`\r\n\x1b[31m[ERROR: ${message}]\x1b[0m\r\n`); // Red color for errors
        },
        onDone: (exitCode) => {
            terminal.current?.writeln(`\r\n\x1b[32m[Script completed with exit code: ${exitCode}]\x1b[0m\r\n`); // Green color
        },
        onConnect: () => {
            terminal.current?.writeln('\x1b[32m[WebSocket Connected]\x1b[0m\r\n');
        },
        onClose: () => {
            terminal.current?.writeln('\r\n\x1b[33m[WebSocket Disconnected]\x1b[0m\r\n');
        },
    });

    // Sync refs with state
    useEffect(() => {
        isConnectedRef.current = isConnected;
        isRunningRef.current = isRunning;
    }, [isConnected, isRunning]);

    // Initialize terminal only once
    useEffect(() => {
        if (!terminalRef.current) return;

        terminal.current = new Terminal({
            cursorBlink: true,
            fontSize: 14,
            fontFamily: 'Menlo, Monaco, "Courier New", monospace',
            theme: {
                background: '#1e1e1e',
                foreground: '#d4d4d4',
                cursor: '#00ff00',
                black: '#000000',
                brightBlack: '#666666',
                red: '#cd3131',
                brightRed: '#f14c4c',
                green: '#0dbc79',
                brightGreen: '#23d18b',
                yellow: '#e5e510',
                brightYellow: '#f5f543',
                blue: '#2472c8',
                brightBlue: '#3b8eea',
                magenta: '#bc3fbc',
                brightMagenta: '#d670d6',
                cyan: '#11a8cd',
                brightCyan: '#29b8db',
                white: '#e5e5e5',
                brightWhite: '#e5e5e5',
            },
        });

        const fitAddon = new FitAddon();
        terminal.current.loadAddon(fitAddon);
        terminal.current.open(terminalRef.current);
        terminal.current.focus();
        fitAddon.fit();

        terminal.current.writeln('Interactive Script Runner Test');
        terminal.current.writeln('Press the "Run Script" button to start\r\n');

        terminal.current.onData((data) => {
            if (isConnectedRef.current && isRunningRef.current) {
                sendInput(data);
            }
        });

        return () => {
            terminal.current?.dispose();
        };
    }, [sendInput]);

    const handleRunScript = () => {
        setHasStarted(true);
        run({
            node_id: '002-cell-stats',
            linux_user: 'john',
            private_key: `-----BEGIN OPENSSH PRIVATE KEY-----
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
-----END OPENSSH PRIVATE KEY-----`,
            environment: 'python-scimap',
            script: 'compute-stats.py',
            args: [
                '--input-selector',
                '../001-raw-input/output/*/*.csv',
                '--output',
                '/nodes/002-cell-stats/output/cell-counts.csv',
            ],
        });
    };

    return (
        <div className="w-full h-screen bg-[#1e1e1e] flex flex-col">
            <div className="p-4 border-b border-gray-700 flex items-center gap-4">
                <h1 className="text-white text-xl font-semibold">Interactive Script Runner</h1>
                <button
                    onClick={handleRunScript}
                    disabled={isRunning}
                    className={`px-4 py-2 rounded ${
                        isRunning
                            ? 'bg-gray-600 cursor-not-allowed'
                            : 'bg-blue-600 hover:bg-blue-700'
                    } text-white transition-colors`}
                >
                    {isRunning ? 'Running...' : 'Run Script'}
                </button>
                {isConnected && (
                    <span className="text-green-400 text-sm">Connected</span>
                )}
                {!isConnected && hasStarted && (
                    <span className="text-red-400 text-sm">Disconnected</span>
                )}
            </div>
            <div className="flex-1 p-4">
                <div ref={terminalRef} className="w-full h-full" />
            </div>
        </div>
    );
}
