"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useEffect, useRef, useState } from "react";

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
        const outputFilename = `data/${crypto.randomUUID()}.png`;
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

            const outputFilename = `data/${crypto.randomUUID()}.csv`;
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
        <div>
            <h1>Proximity Analyze</h1>
            <div>
                <label>Group 1 name: <input type="text" value={group1Name} onChange={e => setGroup1Name(e.target.value)} /></label>
                <input type="file" multiple accept=".csv" onChange={handleGroup1FilesChange} />
            </div>
            <div>
                <label>Group 2 name: <input type="text" value={group2Name} onChange={e => setGroup2Name(e.target.value)} /></label>
                <input type="file" multiple accept=".csv" onChange={handleGroup2FilesChange} />
            </div>

            <button onClick={handleRun}>Run</button>

            <button onClick={handlePlot} disabled={selectedRows.size === 0}>Plot</button>

            <Console ref={consoleRef} />
            
            <table>
                <thead>
                    <tr>
                        {DATA_COLS.map(col => <th key={col}>{colLabel(col)}</th>)}
                        <th>Select</th>
                    </tr>
                </thead>
                <tbody>
                    {tableRows.map((row, i) => (
                        <tr key={i}>
                            {DATA_COLS.map(col => <td key={col}>{row[col]}</td>)}
                            <td>
                                <input type="checkbox" checked={selectedRows.has(i)} onChange={() => toggleRow(i)} />
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>

            {(isPlotting || plotImageUrl) && (
                <div style={{ position: "relative", display: "inline-block" }}>
                    {plotImageUrl && <img src={plotImageUrl} alt="Plot" style={{ opacity: isPlotting ? 0.3 : 1 }} />}
                    {isPlotting && (
                        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center" }}>
                            Generating plot...
                        </div>
                    )}
                </div>
            )}

        </div>
    );
}
