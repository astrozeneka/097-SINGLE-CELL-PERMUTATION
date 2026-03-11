"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useRef, useState } from "react";
import Link from "next/link";

export default function ProximityComputePage() {
    const [headers, setHeaders] = useState<string[][]>([]);
    const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
    const [headerSet, setHeaderSet] = useState<Set<string>>(new Set());
    const [centroidX, setCentroidX] = useState("");
    const [centroidY, setCentroidY] = useState("");
    const [parentArea, setParentArea] = useState("");
    const [parentRegion, setParentRegion] = useState("");
    const [cellType, setCellType] = useState("");
    const [pixelSize, setPixelSize] = useState("0.5");
    const consoleRef = useRef<ConsoleHandle>(null);
    const eventSourceRef = useRef<EventSource | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [downloadFilename, setDownloadFilename] = useState<string | null>(null);

    const handleFileChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
        const files = Array.from(event.target.files ?? []);
        setSelectedFiles(files);
        const results = await Promise.all(
            files.map(file => file.text().then(text => text.split(/\r?\n/)[0].split(",")))
        );
        setHeaders(results);
        // Compute intersection of all header arrays
        const intersection = results.reduce((acc, cur) => acc.filter(h => cur.includes(h)));
        setHeaderSet(new Set(intersection));
    };

    const columnSelect = (value: string, onChange: (v: string) => void, label: string) => (
        <div>
            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">{label}</label>
            <select value={value} onChange={e => onChange(e.target.value)}
                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors">
                <option value="">— select —</option>
                {[...headerSet].map(h => <option key={h} value={h}>{h}</option>)}
            </select>
        </div>
    );

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
            });

        const uploadedFilenames: string[] = [];
        for (const file of selectedFiles) {
            const idx = consoleRef.current!.pushLog(`${file.name} ${bar(0)}`);
            const uploadResult = await uploadFile(file, idx);
            if (!uploadResult.success) {
                consoleRef.current!.pushLog(`ERROR: Failed to upload ${file.name}`, idx);
                setIsProcessing(false);
                return;
            }
            uploadedFilenames.push(uploadResult.path);
        }

        // pass relevant parameters
        const uuid = crypto.randomUUID?.() ?? Math.random().toString(36).slice(2) + Date.now().toString(36);
        const outputFilename = `data/${uuid}.zip`;
        const params = new URLSearchParams({
            script: "proximity-analysis/001_compute.py",
            inputs: uploadedFilenames.join(","),
            centroid_x_col: centroidX,
            centroid_y_col: centroidY,
            parent_area_col: parentArea,
            parent_region_col: parentRegion,
            cell_type_col: cellType,
            pixel_size: pixelSize,
            output: outputFilename,
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
                if (data.exitCode === 0 && data.outputFilename) setDownloadFilename(data.outputFilename);
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
    };

    return (
        <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
            <div className="max-w-6xl mx-auto">
                <header className="mb-6">
                    <Link href="/" className="inline-flex items-center gap-2 text-xs text-slate-400 hover:text-slate-300 transition-colors mb-4">
                        ← Back to home
                    </Link>
                    <div className="flex items-center gap-4">
                        <h1 className="text-2xl font-light text-slate-100 tracking-wide">
                            Proximity Analysis (Step 1): Compute
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
                    <div>
                        <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Input Files</label>
                        <input type="file" multiple accept=".csv" onChange={handleFileChange}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700" />
                    </div>
                    <div className="grid grid-cols-2 gap-4">
                        {columnSelect(centroidX, setCentroidX, "Centroid X Column")}
                        {columnSelect(centroidY, setCentroidY, "Centroid Y Column")}
                        {columnSelect(parentRegion, setParentRegion, "Parent Region Column")}
                        {columnSelect(parentArea, setParentArea, "Parent Area Column")}
                        {columnSelect(cellType, setCellType, "Cell Type Column")}
                        <div>
                            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Pixel Size (µm)</label>
                            <input type="number" value={pixelSize} onChange={e => setPixelSize(e.target.value)} step="0.1" min="0"
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors" />
                        </div>
                    </div>
                    <div className="flex gap-3 items-center">
                        <button onClick={handleRun} disabled={isProcessing}
                            className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors">
                            {isProcessing ? "Running..." : "Run"}
                        </button>
                        {downloadFilename && !isProcessing && (
                            <a href={`/api/download-file?filename=${downloadFilename}`} download
                                className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors">
                                ↓ Download ZIP
                            </a>
                        )}
                    </div>
                </div>

                <div className="grid grid-cols-2 gap-6">
                    <div>
                        <Console ref={consoleRef} />
                    </div>
                    <div>
                        <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider mb-3">Output File</h2>
                        <div className="bg-slate-900/50 backdrop-blur border-l border-slate-800 h-[calc(100vh-480px)] overflow-auto flex items-center justify-center p-3">
                            {!isProcessing && downloadFilename ? (
                                <div className="text-slate-300 text-center">
                                    <div className="text-emerald-400 text-lg mb-2">✓</div>
                                    <div className="text-sm">Output ready for download</div>
                                    <div className="text-xs text-slate-500 mt-1">{downloadFilename.replace(/^data\//, "")}</div>
                                </div>
                            ) : (
                                <div className="text-slate-600 italic text-sm">No output generated yet...</div>
                            )}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}
