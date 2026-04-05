'use client';

import { FitAddon } from "@xterm/addon-fit";
import { Terminal } from "@xterm/xterm";
import { useEffect, useRef, forwardRef, useImperativeHandle } from "react";
import { useScriptRunner, ScriptRunnerConfig } from "@/hooks/use-script-runner";
import '@xterm/xterm/css/xterm.css';

export interface RunnerTerminalRef {
    run: (config: ScriptRunnerConfig) => void;
    isRunning: boolean;
    isConnected: boolean;
}

const RunnerTerminal = forwardRef<RunnerTerminalRef>((_props, ref) => {
    const terminalRef = useRef<HTMLDivElement>(null);
    const terminal = useRef<Terminal | null>(null);
    const isConnectedRef = useRef(false);
    const isRunningRef = useRef(false);

    const { run, sendInput, isConnected, isRunning } = useScriptRunner({
        onStdout: (data) => {
            terminal.current?.write(data);
        },
        onStderr: (data) => {
            terminal.current?.write(`\x1b[31m${data}\x1b[0m`);
        },
        onStatus: (message) => {
            terminal.current?.writeln(`\r\n\x1b[36m[${message}]\x1b[0m\r\n`);
        },
        onError: (message) => {
            terminal.current?.writeln(`\r\n\x1b[31m[ERROR: ${message}]\x1b[0m\r\n`);
        },
        onDone: (exitCode) => {
            terminal.current?.writeln(`\r\n\x1b[32m[Script completed with exit code: ${exitCode}]\x1b[0m\r\n`);
        },
        onConnect: () => {
            terminal.current?.writeln('\x1b[32m[WebSocket Connected]\x1b[0m\r\n');
        },
        onClose: () => {
            terminal.current?.writeln('\r\n\x1b[33m[WebSocket Disconnected]\x1b[0m\r\n');
        },
    });

    useEffect(() => {
        isConnectedRef.current = isConnected;
        isRunningRef.current = isRunning;
    }, [isConnected, isRunning]);

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

        terminal.current.onData((data) => {
            if (isConnectedRef.current && isRunningRef.current) {
                sendInput(data);
            }
        });

        return () => {
            terminal.current?.dispose();
        };
    }, [sendInput]);

    useImperativeHandle(ref, () => ({
        run,
        isRunning,
        isConnected
    }));

    return (
        <div className="w-full h-screen bg-[#1e1e1e] p-4">
            <div ref={terminalRef} className="w-full h-full" />
        </div>
    );
});

RunnerTerminal.displayName = 'RunnerTerminal';

export default RunnerTerminal;