"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useEffect, useRef, useState } from "react";
import Link from "next/link";

const DATA_COLS = ["Phenotype A", "Phenotype B", "Mean Group 1", "Mean Group 2", "Mann-Whitney U", "p-value"];

export default function ProximityComputePage() {

    const [group1Name, setGroup1Name] = useState("Group 1");
    const [group2Name, setGroup2Name] = useState("Group 2");
    const [group1Files, setGroup1Files] = useState<File[]>([]);
    const [group2Files, setGroup2Files] = useState<File[]>([]);
    const consoleRef = useRef<ConsoleHandle>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const eventSourceRef = useRef<EventSource | null>(null);
    const [downloadCSVFilename, setDownloadCSVFilename] = useState<string | null>(null);
    const [tableRows, setTableRows] = useState<Record<string, string>[]>([]);
    const [selectedRows, setSelectedRows] = useState<Set<number>>(new Set());
    const [plotImageUrl, setPlotImageUrl] = useState<string | null>(null);
    const [isPlotting, setIsPlotting] = useState(false);

    const colLabel = (col: string) => {
        if (col === "Mean Group 1") return `Mean ${group1Name}`;
        if (col === "Mean Group 2") return `Mean ${group2Name}`;
        return col;
    };

    const handleFilesChange = (event: React.ChangeEvent<HTMLInputElement>, mutator: React.Dispatch<React.SetStateAction<File[]>>) => {
        const files = Array.from(event.target.files ?? []);
        mutator(files);
    };

    const handleGroup1FilesChange = (e: React.ChangeEvent<HTMLInputElement>) => handleFilesChange(e, setGroup1Files);
    const handleGroup2FilesChange = (e: React.ChangeEvent<HTMLInputElement>) => handleFilesChange(e, setGroup2Files);

    useEffect(() => {
        if (!downloadCSVFilename) return;
        const filename = downloadCSVFilename.replace(/^data\//, "");
        fetch(`/api/download-file?filename=${filename}`)
            .then(r => r.text())
            .then(text => {
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
                const lines = text.trim().split(/\r?\n/);
                const headers = parseCSVLine(lines[0]);
                const rows = lines.slice(1).map(line => {
                    const values = parseCSVLine(line);
                    return Object.fromEntries(headers.map((h, i) => [h, values[i] ?? ""]));
                });
                setTableRows(rows);
                setSelectedRows(new Set());
            });
    }, [downloadCSVFilename]);

    const toggleRow = (i: number) => {
        setSelectedRows(prev => {
            const next = new Set(prev);
            next.has(i) ? next.delete(i) : next.add(i);
            return next;
        });
    };

    const handlePlot = () => {
        if (!downloadCSVFilename || selectedRows.size === 0) return;
        setIsPlotting(true);
        const tuples = [...selectedRows]
            .map(i => `${tableRows[i]['Phenotype A']}|${tableRows[i]['Phenotype B']}`)
            .join(",");
        const uuid = crypto.randomUUID?.() ?? Math.random().toString(36).slice(2) + Date.now().toString(36);
        const outputFilename = `data/${uuid}.png`;
        const params = new URLSearchParams({
            script: "proximity-analysis/003_plot.py",
            input: downloadCSVFilename,
            tuples,
            group1_name: group1Name,
            group2_name: group2Name,
            output: outputFilename,
        });
        const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
        eventSource.onmessage = (event) => {
            const data = JSON.parse(event.data);
            if (data.type === "complete" && data.exitCode === 0) {
                const filename = outputFilename.replace(/^data\//, "");
                fetch(`/api/download-file?filename=${filename}`)
                    .then(r => r.blob())
                    .then(blob => { setPlotImageUrl(URL.createObjectURL(blob)); setIsPlotting(false); });
            }
            if (data.type === "complete" || data.type === "error") { eventSource.close(); setIsPlotting(false); }
        };
    };

    const handleRun = async () => {
        setIsProcessing(true);
        consoleRef.current?.clearLogs();

        const bar = (pct: number) => {
            const filled = Math.round(pct / 5);
            return `[${'#'.repeat(filled)}${'-'.repeat(20 - filled)}] ${pct}%`;
        };
        const uploadFile = (file: File, idx: number): Promise<any> =>
            new Promise((resolve, reject) => {
                const xhr = new XMLHttpRequest();
                const formData = new FormData();
                formData.append("file", file);
                xhr.upload.onprogress = (e) => {
                    if (e.lengthComputable)
                        consoleRef.current!.pushLog(`${file.name} ${bar(Math.round(e.loaded / e.total * 100))}`, idx);
                };
                xhr.onload = () => resolve(JSON.parse(xhr.responseText));
                xhr.onerror = () => reject(new Error("Upload failed"));
                xhr.open("POST", "/api/upload-file");
                xhr.send(formData);
            })

            const uploadedGroup1FilePaths: string[] = [];
            const uploadedGroup2FilePaths: string[] = [];

            const uploadedFilePaths: string[] = [];
            for (const [files, uploadedPaths] of [[group1Files, uploadedGroup1FilePaths], [group2Files, uploadedGroup2FilePaths]] as const) {
                for (const file of files) {
                    const idx = consoleRef.current!.pushLog(`${file.name} ${bar(0)}`);
                    const uploadResult = await uploadFile(file, idx);
                    if (!uploadResult.success) {
                        consoleRef.current!.pushLog(`ERROR: Failed to upload ${file.name}`, idx);
                        setIsProcessing(false);
                        return;
                    }
                    uploadedPaths.push(uploadResult.path);
                }
            }

            const uuid = crypto.randomUUID?.() ?? Math.random().toString(36).slice(2) + Date.now().toString(36);
            const outputFilename = `data/${uuid}.csv`;
            const params = new URLSearchParams({
                script: "proximity-analysis/002_analyze.py",
                group1_inputs: uploadedGroup1FilePaths.join(","),
                group2_inputs: uploadedGroup2FilePaths.join(","),
                output: outputFilename
            });
            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
            eventSourceRef.current = eventSource;
            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);
                if (data.type === "acknowledgment") {
                    consoleRef.current?.pushLog(data.message);
                } else if (data.type === "stdout") {
                    consoleRef.current?.pushLog(data.content);
                } else if (data.type === "stderr") {
                    consoleRef.current?.pushLog(`ERROR: ${data.content}`);
                } else if (data.type === "complete") {
                    consoleRef.current?.pushLog(`Process completed with exit code: ${data.exitCode}`);
                    if (data.exitCode === 0 && data.outputFilename) setDownloadCSVFilename(data.outputFilename);
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    consoleRef.current?.pushLog(`ERROR: ${data.message}`);
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                }
            };
            eventSource.onerror = () => {
                consoleRef.current?.pushLog("Connection error");
                eventSource.close();
                eventSourceRef.current = null;
                setIsProcessing(false);
            };
    }

    return (
        <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
            <div className="max-w-6xl mx-auto">
                <header className="mb-6">
                    <Link href="/" className="inline-flex items-center gap-2 text-xs text-slate-400 hover:text-slate-300 transition-colors mb-4">
                        ← Back to home
                    </Link>
                    <div className="flex items-center gap-4">
                        <h1 className="text-2xl font-light text-slate-100 tracking-wide">
                            Proximity Analysis (Step 2): Analyze
                        </h1>
                        {isProcessing && (
                            <span className="text-xs text-red-400 font-medium uppercase tracking-wider animate-pulse">
                                Do not close this page
                            </span>
                        )}
                    </div>
                    <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
                </header>

                <div className="mb-6 space-y-4">
                    <div className="grid grid-cols-2 gap-4">
                        <div>
                            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Group 1 Name</label>
                            <input type="text" value={group1Name} onChange={e => setGroup1Name(e.target.value)}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors" />
                        </div>
                        <div>
                            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Group 1 Files</label>
                            <input type="file" multiple accept=".csv" onChange={handleGroup1FilesChange}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700" />
                        </div>
                    </div>
                    <div className="grid grid-cols-2 gap-4">
                        <div>
                            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Group 2 Name</label>
                            <input type="text" value={group2Name} onChange={e => setGroup2Name(e.target.value)}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors" />
                        </div>
                        <div>
                            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Group 2 Files</label>
                            <input type="file" multiple accept=".csv" onChange={handleGroup2FilesChange}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700" />
                        </div>
                    </div>
                    <div className="flex gap-3 items-center">
                        <button onClick={handleRun} disabled={isProcessing}
                            className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors">
                            {isProcessing ? "Running..." : "Run"}
                        </button>
                        <button onClick={handlePlot} disabled={selectedRows.size === 0 || isPlotting}
                            className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors">
                            {isPlotting ? "Plotting..." : "Plot"}
                        </button>
                        {downloadCSVFilename && !isProcessing && (
                            <a href={`/api/download-file?filename=${downloadCSVFilename.replace(/^data\//, "")}`} download
                                className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors">
                                ↓ Download CSV
                            </a>
                        )}
                    </div>
                </div>

                <div className="grid grid-cols-2 gap-6 mb-6">
                    <div>
                        <Console ref={consoleRef} />
                    </div>
                    <div>
                        <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider mb-3">Results</h2>
                        <div className="bg-slate-900/50 border-l border-slate-800 overflow-y-auto overflow-x-hidden h-[300px]">
                            {tableRows.length > 0 ? (
                                <table className="w-full table-fixed text-xs">
                                    <thead className="sticky top-0 bg-slate-900">
                                        <tr>
                                            {DATA_COLS.map(col => (
                                                <th key={col} className="px-1.5 py-1.5 text-left text-slate-400 font-medium uppercase tracking-wider truncate">
                                                    {colLabel(col)}
                                                </th>
                                            ))}
                                            <th className="px-1.5 py-1.5 text-center text-slate-400 font-medium w-7">✓</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        {tableRows.map((row, i) => (
                                            <tr key={i} onClick={() => toggleRow(i)}
                                                className={`cursor-pointer border-t border-slate-800 hover:bg-slate-800/50 ${selectedRows.has(i) ? "bg-slate-800/70" : ""}`}>
                                                {DATA_COLS.map(col => (
                                                    <td key={col} className="px-1.5 py-1 text-slate-300 truncate">{row[col]}</td>
                                                ))}
                                                <td className="px-1.5 py-1 text-center">
                                                    <input type="checkbox" checked={selectedRows.has(i)} onChange={() => toggleRow(i)} onClick={e => e.stopPropagation()} className="accent-teal-400" />
                                                </td>
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            ) : (
                                <div className="h-full flex items-center justify-center text-slate-600 italic text-sm">No results yet...</div>
                            )}
                        </div>
                    </div>
                </div>

                {(isPlotting || plotImageUrl) && (
                    <div>
                        <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider mb-3">Figure</h2>
                        <div className="bg-slate-900/50 border-l border-slate-800 p-4 relative">
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
    );
}
