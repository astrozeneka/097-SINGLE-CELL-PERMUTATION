import { FitAddon } from "@xterm/addon-fit";
import { Terminal } from "@xterm/xterm";
import { useEffect, useRef } from "react";


export default function RunnerTerminal() {
    const terminalRef = useRef<HTMLDivElement>(null);
    const terminal = useRef<Terminal | null>(null);
    const ws = useRef<WebSocket | null>(null);

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
        fitAddon.fit();

    }, []);

    const start = (linux_user: string, private_key: string, environment: string, script: string, args?: string[]) => {

    };


    return (
        <div className="w-full h-screen bg-[#1e1e1e] p-4">
            <div ref={terminalRef} className="w-full h-full" />
        </div>
    );
}