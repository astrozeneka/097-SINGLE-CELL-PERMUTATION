"use client";

import { forwardRef, useImperativeHandle, useState } from "react";

export interface LogDisplayHandle {
    pushLog: (message: string) => void;
    clearLogs: () => void;
}

export const LogDisplay = forwardRef<LogDisplayHandle>((_, ref) => {
    const [logs, setLogs] = useState<string[]>([]);

    useImperativeHandle(ref, () => ({
        pushLog: (message: string) => {
            setLogs((prev) => [...prev, message]);
        },
        clearLogs: () => {
            setLogs([]);
        },
    }));

    return (
        <div className="bg-slate-900/50 backdrop-blur font-mono text-xs p-3 h-[calc(100vh-480px)] overflow-y-auto flex flex-col-reverse border-l border-slate-800">
            {logs.length === 0 ? (
                <div className="text-slate-600 italic">Waiting for input...</div>
            ) : (
                [...logs].reverse().map((log, index) => (
                    <div key={index} className="text-slate-300 leading-relaxed py-0.5">
                        {log}
                    </div>
                ))
            )}
        </div>
    );
});

LogDisplay.displayName = "LogDisplay";
