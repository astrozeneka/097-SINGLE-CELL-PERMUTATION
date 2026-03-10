"use client";

import { forwardRef, useImperativeHandle, useState } from "react";

export interface ConsoleHandle {
    pushLog: (message: string) => void;
    clearLogs: () => void;
}

const Console = forwardRef<ConsoleHandle>(function Console(_, ref) {
    const [logs, setLogs] = useState<string[]>([]);

    useImperativeHandle(ref, () => ({
        pushLog: (message: string) => setLogs(prev => [...prev, message]),
        clearLogs: () => setLogs([]),
    }));

    return (
        <div>
            <div className="flex justify-between items-center mb-3">
                <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
                    Console Output
                </h2>
            </div>
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
        </div>
    );
});

export default Console;
