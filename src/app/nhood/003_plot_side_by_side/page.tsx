"use client";

import { FormEvent, useEffect, useRef, useState } from "react";
import Link from "next/link";
import { LogDisplay, LogDisplayHandle } from "../components/log-display";


export default function NHoodPlotSideBySide() {
    const [file, setFile] = useState<File | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [outputFilename, setOutputFilename] = useState<string | null>(null);
    const logRef = useRef<LogDisplayHandle>(null);
    const eventSourceRef = useRef<EventSource | null>(null);

    const handleSubmit = async (e: FormEvent) => {
        e.preventDefault();
        if (!file) {
            alert("Please select a file");
            return;
        }
        setIsProcessing(true);
        logRef.current?.clearLogs();
        setOutputFilename(null);

        try {
            // Upload file
            const formData = new FormData();
            formData.append("file", file);

            logRef.current?.pushLog("Uploading file...");
            const uploadResponse = await fetch("/api/upload-file", {
                method: "POST",
                body: formData,
            });

            const uploadResult = await uploadResponse.json();
            if (!uploadResult.success) {
                setIsProcessing(false);
                logRef.current?.pushLog("File upload failed");
                throw new Error("File upload failed");
            }

            logRef.current?.pushLog("File uploaded successfully");

            // Run the plot script
            logRef.current?.pushLog("Starting script...");
            const outputName = uploadResult.filename.replace('.csv', '_plot.png');
            const params = new URLSearchParams({
                input: `data/${uploadResult.filename}`,
                output: `data/${outputName}`,
                script: "nhood/003_plot_side_by_side.py",
            });
            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
            eventSourceRef.current = eventSource;

            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);

                if (data.type === "acknowledgment") {
                    logRef.current?.pushLog(data.message);
                } else if (data.type === "stdout") {
                    logRef.current?.pushLog(data.content);
                } else if (data.type === "stderr") {
                    logRef.current?.pushLog(`ERROR: ${data.content}`);
                } else if (data.type === "complete") {
                    logRef.current?.pushLog(`Process completed with exit code: ${data.exitCode}`);
                    if (data.exitCode === 0 && data.outputFilename) {
                        setOutputFilename(data.outputFilename.replace('data/', ''));
                    }
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    logRef.current?.pushLog(`ERROR: ${data.message}`);
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                }
            };
            eventSource.onerror = () => {
                logRef.current?.pushLog("Connection error");
                eventSource.close();
                eventSourceRef.current = null;
                setIsProcessing(false);
            };

        } catch (error: any) {
            logRef.current?.pushLog(`Error: ${error.message}`);
            setIsProcessing(false);
        }

    }

    useEffect(() => {
        return () => {
            if (eventSourceRef.current) {
                eventSourceRef.current.close();
                eventSourceRef.current = null;
            }
        };
    }, []);

    return (
        <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
            <div className="max-w-5xl mx-auto">
                <header className="mb-6">
                    <Link
                        href="/"
                        className="inline-flex items-center gap-2 text-xs text-slate-400 hover:text-slate-300 transition-colors mb-4"
                    >
                        ← Back to home
                    </Link>
                    <div className="flex items-center gap-4">
                        <h1 className="text-2xl font-light text-slate-100 tracking-wide">
                            Neighborhood Analysis (Step 3): Visualize
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
                        <div className="flex justify-between items-center mb-1.5">
                            <label className="text-xs text-slate-400 uppercase tracking-wider">
                                Input File
                            </label>
                            {/*<a
                                href="/sample-datas/nhood-cluster-input-data.csv"
                                download
                                className="text-xs text-teal-400 hover:text-teal-300 transition-colors"
                            >
                                ↓ Download sample input
                            </a>*/}
                        </div>
                        <input
                            type="file"
                            accept=".csv"
                            onChange={(e) => setFile(e.target.files?.[0] || null)}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                        />
                    </div>

                    <button
                        type="submit"
                        disabled={isProcessing}
                        className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
                    >
                        {isProcessing ? "Running..." : "Run"}
                    </button>
                </form>

                <div className="grid grid-cols-2 gap-6">
                    <div>
                        <div className="flex justify-between items-center mb-3">
                            <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
                                Console Output
                            </h2>
                        </div>
                        <LogDisplay ref={logRef} />
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
                                    ↓ Download PNG
                                </a>
                            )}
                        </div>
                        <div className="bg-slate-900/50 backdrop-blur border-l border-slate-800 h-[calc(100vh-480px)] overflow-auto flex items-center justify-center p-3">
                            {!isProcessing && outputFilename ? (
                                <div className="text-slate-300 text-center">
                                    <div className="text-emerald-400 text-lg mb-2">✓</div>
                                    <div className="text-sm">Plot ready for download</div>
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
