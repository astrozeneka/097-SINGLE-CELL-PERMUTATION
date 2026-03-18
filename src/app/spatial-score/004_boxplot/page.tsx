"use client";

import { useState, useRef } from "react";
import Console, { ConsoleHandle } from "@/app/components/console";
import Link from "next/link";

const parseCSVLine = (line: string) => {
    const fields: string[] = [];
    let cur = "", inQuote = false;
    for (const ch of line) {
        if (ch === '"') { inQuote = !inQuote; }
        else if (ch === ',' && !inQuote) { fields.push(cur); cur = ""; }
        else { cur += ch; }
    }
    fields.push(cur);
    return fields;
};

export default function Boxplot() {
    const [tableRows, setTableRows] = useState<Record<string, string>[]>([]);
    const [displayCols, setDisplayCols] = useState<string[]>([]);
    const [selectedRows, setSelectedRows] = useState<Set<number>>(new Set());
    const [inputFile, setInputFile] = useState<File | null>(null);
    const [serverPath, setServerPath] = useState<string | null>(null);
    const [plotImageUrl, setPlotImageUrl] = useState<string | null>(null);
    const [isPlotting, setIsPlotting] = useState(false);
    const consoleRef = useRef<ConsoleHandle>(null);

    const handleFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        setInputFile(file);
        setServerPath(null);
        setPlotImageUrl(null);
        setSelectedRows(new Set());

        const text = await file.text();
        const lines = text.trim().split(/\r?\n/);
        const headers = parseCSVLine(lines[0]);
        setDisplayCols(headers.filter(h => !h.endsWith('_vals')));
        const rows = lines.slice(1).map(line => {
            const values = parseCSVLine(line);
            return Object.fromEntries(headers.map((h, i) => [h, values[i] ?? ""]));
        });
        setTableRows(rows);

        // Upload immediately so plot is fast
        const formData = new FormData();
        formData.append("file", file);
        const result = await fetch("/api/upload-file", { method: "POST", body: formData }).then(r => r.json());
        if (result.success) setServerPath(result.path);
    };

    const toggleRow = (i: number) => {
        setSelectedRows(prev => {
            const next = new Set(prev);
            next.has(i) ? next.delete(i) : next.add(i);
            return next;
        });
    };

    const handlePlot = () => {
        if (!serverPath || selectedRows.size === 0) return;
        setIsPlotting(true);
        consoleRef.current?.clearLogs();

        const tuples = [...selectedRows]
            .map(i => `${tableRows[i]['target']}|${tableRows[i]['effector']}`)
            .join(",");
        const uuid = Math.random().toString(36).slice(2) + Date.now().toString(36);
        const outputFilename = `${uuid}.png`;

        const params = new URLSearchParams({
            script: "spatial-score/004_boxplot.py",
            input: serverPath,
            tuples,
            output: `data/${outputFilename}`,
        });

        const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
        eventSource.onmessage = (event) => {
            const data = JSON.parse(event.data);
            if (data.type === "stdout") consoleRef.current?.pushLog(data.content);
            if (data.type === "stderr") consoleRef.current?.pushLog(`ERROR: ${data.content}`);
            if (data.type === "complete") {
                consoleRef.current?.pushLog(`Exit code: ${data.exitCode}`);
                if (data.exitCode === 0) {
                    fetch(`/api/download-file?filename=${outputFilename}`)
                        .then(r => r.blob())
                        .then(blob => { setPlotImageUrl(URL.createObjectURL(blob)); setIsPlotting(false); });
                } else {
                    setIsPlotting(false);
                }
                eventSource.close();
            }
            if (data.type === "error") {
                consoleRef.current?.pushLog(`ERROR: ${data.message}`);
                eventSource.close();
                setIsPlotting(false);
            }
        };
        eventSource.onerror = () => { eventSource.close(); setIsPlotting(false); };
    };

    return (
        <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
            <div className="max-w-6xl mx-auto">
                <header className="mb-6">
                    <Link href="/" className="inline-flex items-center gap-2 text-xs text-slate-400 hover:text-slate-300 transition-colors mb-4">
                        ← Back to home
                    </Link>
                    <h1 className="text-2xl font-light text-slate-100 tracking-wide">Spatial Score: Boxplot</h1>
                    <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
                </header>

                <div className="mb-6 space-y-4">
                    <div>
                        <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Input File</label>
                        <input type="file" accept=".csv" onChange={handleFileChange}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700" />
                    </div>

                    {tableRows.length > 0 && (
                        <>
                            <div className="bg-slate-900/50 border-l border-slate-800 overflow-y-auto overflow-x-auto h-[300px]">
                                <table className="w-full table-fixed text-xs">
                                    <thead className="sticky top-0 bg-slate-900">
                                        <tr>
                                            {displayCols.map(col => (
                                                <th key={col} className="px-1.5 py-1.5 text-left text-slate-400 font-medium uppercase tracking-wider truncate">
                                                    {col}
                                                </th>
                                            ))}
                                            <th className="px-1.5 py-1.5 text-center text-slate-400 font-medium w-7">✓</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        {tableRows.map((row, i) => (
                                            <tr key={i} onClick={() => toggleRow(i)}
                                                className={`cursor-pointer border-t border-slate-800 hover:bg-slate-800/50 ${selectedRows.has(i) ? "bg-slate-800/70" : ""}`}>
                                                {displayCols.map(col => (
                                                    <td key={col} className="px-1.5 py-1 text-slate-300 truncate">{row[col]}</td>
                                                ))}
                                                <td className="px-1.5 py-1 text-center">
                                                    <input type="checkbox" checked={selectedRows.has(i)} onChange={() => toggleRow(i)}
                                                        onClick={e => e.stopPropagation()} className="accent-teal-400" />
                                                </td>
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            </div>

                            <div className="flex gap-3 items-center">
                                <button onClick={handlePlot} disabled={selectedRows.size === 0 || isPlotting || !serverPath}
                                    className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors">
                                    {isPlotting ? "Plotting..." : "Plot"}
                                </button>
                            </div>
                        </>
                    )}
                </div>

                <div className="grid grid-cols-2 gap-6 mb-6">
                    <Console ref={consoleRef} />
                    {(isPlotting || plotImageUrl) && (
                        <div>
                            <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider mb-3">Figure</h2>
                            <div className="bg-white border-l border-slate-800 p-4 relative">
                                {plotImageUrl && <img src={plotImageUrl} alt="Plot" className={`max-w-full ${isPlotting ? "opacity-30" : ""}`} />}
                                {isPlotting && (
                                    <div className="absolute inset-0 flex items-center justify-center text-slate-400 text-sm">
                                        Generating plot...
                                    </div>
                                )}
                            </div>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}
