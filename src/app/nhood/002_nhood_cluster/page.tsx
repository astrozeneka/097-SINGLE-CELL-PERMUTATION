"use client";

import { FormEvent, useRef, useState } from "react";
import { LogDisplay, LogDisplayHandle } from "../components/log-display";


export default function NHoodClusterPage() {
    const [file, setFile] = useState<File | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [outputFilename, setOutputFilename] = useState<string | null>(null);
    const logRef = useRef<LogDisplayHandle>(null);

    // The parameters for the submission
    const [nMotifs, setNMotifs] = useState(20);
    const [knnSeed, setKnnSeed] = useState(42);

    // 
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

            // Run the clustering process
            // Begin by preparing the request
            logRef.current?.pushLog("Starting script...");
            const outputName = uploadResult.filename.replace('.csv', '_withmotifs.csv');
            const params = new URLSearchParams({
                input: `data/${uploadResult.filename}`,
                output: `data/${outputName}`,
                n_motifs: nMotifs.toString(),
                knn_seed: knnSeed.toString(),
            });
            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
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
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    logRef.current?.pushLog(`ERROR: ${data.message}`);
                    eventSource.close();
                    setIsProcessing(false);
                }
            };
            eventSource.onerror = () => {
                logRef.current?.pushLog("Connection error");
                eventSource.close();
                setIsProcessing(false);
            };

        } catch (error: any) {
            logRef.current?.pushLog(`Error: ${error.message}`);
            setIsProcessing(false);
        } finally {
            setIsProcessing(false);
        }

    }

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
                            N Motifs
                        </label>
                        <input
                            type="number"
                            value={nMotifs}
                            onChange={(e) => setNMotifs(parseInt(e.target.value))}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
                        />
                    </div>

                    <div>
                        <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                            KNN Seed
                        </label>
                        <input
                            type="number"
                            value={knnSeed}
                            onChange={(e) => setKnnSeed(parseInt(e.target.value))}
                            className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
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

                <div className="mb-6">
                    <div className="flex justify-between items-center mb-3">
                        <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
                            Console Output
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
                    <LogDisplay ref={logRef} />
                </div>

            </div>
        </div>

    )
}