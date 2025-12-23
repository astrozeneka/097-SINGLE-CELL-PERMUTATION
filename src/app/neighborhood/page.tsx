"use client";

import { useState, FormEvent } from "react";

export default function NeighborhoodPage() {
    const [file, setFile] = useState<File | null>(null);
    const [logs, setLogs] = useState<string[]>([]);
    const [isProcessing, setIsProcessing] = useState(false);
    const [outputFilename, setOutputFilename] = useState<string | null>(null);
    const [method, setMethod] = useState<string>("knn");
    const [knnCount, setKnnCount] = useState<number>(20);
    const [radius, setRadius] = useState<number>(50.0);

    const handleSubmit = async (e: FormEvent) => {
        e.preventDefault();
        if (!file) {
            alert("Please select a file");
            return;
        }

        setIsProcessing(true);
        setLogs([]);
        setOutputFilename(null);

        try {
            // Upload file
            const formData = new FormData();
            formData.append("file", file);

            setLogs((prev) => [...prev, "Uploading file..."]);
            const uploadResponse = await fetch("/api/upload-file", {
                method: "POST",
                body: formData,
            });

            const uploadResult = await uploadResponse.json();
            if (!uploadResult.success) {
                throw new Error("File upload failed");
            }

            setLogs((prev) => [...prev, `File uploaded: ${uploadResult.filename}`]);

            // Run neighborhood analysis with SSE
            setLogs((prev) => [...prev, "Starting neighborhood analysis..."]);
            const outputName = uploadResult.filename.replace('.csv', '_cnhood.csv');
            const params = new URLSearchParams({
                input: `data/${uploadResult.filename}`,
                output: `data/${outputName}`,
                script: "nhood/001_neighborhood.py",
                method: method,
                ...(method === "knn" && { "knn-count": knnCount.toString() }),
                ...(method === "radius" && { radius: radius.toString() })
            });
            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);

            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);

                if (data.type === "acknowledgment") {
                    setLogs((prev) => [...prev, data.message]);
                } else if (data.type === "stdout") {
                    setLogs((prev) => [...prev, data.content]);
                } else if (data.type === "stderr") {
                    setLogs((prev) => [...prev, `ERROR: ${data.content}`]);
                } else if (data.type === "complete") {
                    setLogs((prev) => [...prev, `Process completed with exit code: ${data.exitCode}`]);
                    if (data.exitCode === 0 && data.outputFilename) {
                        setOutputFilename(data.outputFilename.replace('data/', ''));
                    }
                    eventSource.close();
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    setLogs((prev) => [...prev, `ERROR: ${data.message}`]);
                    eventSource.close();
                    setIsProcessing(false);
                }
            };

            eventSource.onerror = () => {
                setLogs((prev) => [...prev, "Connection error"]);
                eventSource.close();
                setIsProcessing(false);
            };
        } catch (error: any) {
            setLogs((prev) => [...prev, `Error: ${error.message}`]);
            setIsProcessing(false);
        }
    };

    return (
        <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
            <div className="max-w-5xl mx-auto">
                <header className="mb-6">
                    <div className="flex items-center gap-4">
                        <h1 className="text-2xl font-light text-slate-100 tracking-wide">
                            Neighborhood Analysis
                        </h1>
                        {isProcessing && (
                            <span className="text-xs text-red-400 font-medium uppercase tracking-wider animate-pulse">
                                Do not close this page
                            </span>
                        )}
                    </div>
                    <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
                </header>

                <form onSubmit={handleSubmit} className="mb-6 space-y-4">
                    <div>
                        <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                            Input File
                        </label>
                        <input
                            type="file"
                            accept=".csv"
                            onChange={(e) => setFile(e.target.files?.[0] || null)}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                        />
                    </div>

                    <div>
                        <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                            Method
                        </label>
                        <select
                            value={method}
                            onChange={(e) => setMethod(e.target.value)}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
                        >
                            <option value="knn">KNN</option>
                            <option value="radius">Radius</option>
                        </select>
                    </div>

                    {method === "knn" && (
                        <div>
                            <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                                KNN Count
                            </label>
                            <input
                                type="number"
                                value={knnCount}
                                onChange={(e) => setKnnCount(parseInt(e.target.value))}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
                            />
                        </div>
                    )}

                    {method === "radius" && (
                        <div>
                            <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                                Radius (px) - <span className="text-red-500">make sure to properly convert micrometers to pixels if needed</span>
                            </label>
                            <input
                                type="number"
                                step="5"
                                value={radius}
                                onChange={(e) => setRadius(parseFloat(e.target.value))}
                                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
                            />
                        </div>
                    )}

                    <button
                        type="submit"
                        disabled={isProcessing}
                        className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
                    >
                        {isProcessing ? "Running..." : "Generate"}
                    </button>
                </form>

                <div className="grid grid-cols-2 gap-6">
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

                    <div>
                        <div className="flex justify-between items-center mb-3">
                            <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
                                Output File
                            </h2>
                            {!isProcessing && outputFilename && (
                                <a
                                    href={`/api/download-file?filename=${outputFilename}`}
                                    download
                                    className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors"
                                >
                                    ↓ Download CSV
                                </a>
                            )}
                        </div>
                        <div className="bg-slate-900/50 backdrop-blur border-l border-slate-800 h-[calc(100vh-480px)] overflow-auto flex items-center justify-center p-3">
                            {!isProcessing && outputFilename ? (
                                <div className="text-slate-300 text-center">
                                    <div className="text-emerald-400 text-lg mb-2">✓</div>
                                    <div className="text-sm">File ready for download</div>
                                    <div className="text-xs text-slate-500 mt-1">{outputFilename}</div>
                                </div>
                            ) : (
                                <div className="text-slate-600 italic">No file generated yet...</div>
                            )}
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}