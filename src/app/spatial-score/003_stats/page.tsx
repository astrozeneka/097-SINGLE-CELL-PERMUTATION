"use client";

import { useState, useRef } from "react";
import Console, { ConsoleHandle } from "@/app/components/console";
import Link from "next/link";

const MIN_GROUPS = 2;
const MAX_GROUPS = 4;

type Group = { name: string; files: File[] };

const emptyGroup = (): Group => ({ name: "", files: [] });

export default function Boxplot() {
    const [groups, setGroups] = useState<Group[]>([emptyGroup(), emptyGroup()]);
    const [isProcessing, setIsProcessing] = useState(false);
    const [downloadFilename, setDownloadFilename] = useState<string | null>(null);
    const consoleRef = useRef<ConsoleHandle>(null);
    const eventSourceRef = useRef<EventSource | null>(null);

    const addGroup = () => setGroups(g => [...g, emptyGroup()]);
    const removeGroup = (i: number) => setGroups(g => g.filter((_, idx) => idx !== i));

    const updateName = (i: number, name: string) =>
        setGroups(g => g.map((grp, idx) => idx === i ? { ...grp, name } : grp));

    const updateFiles = (i: number, files: FileList | null) =>
        setGroups(g => g.map((grp, idx) => idx === i ? { ...grp, files: files ? Array.from(files) : [] } : grp));

    const runValid = groups.every(g => g.files.length > 0);

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

    const handleRun = async () => {
        setIsProcessing(true);
        setDownloadFilename(null);
        consoleRef.current?.clearLogs();

        const groupUploadedPaths: Record<string, string[]> = {};
        for (let i = 0; i < groups.length; i++) {
            const groupKey = `group${i + 1}`;
            groupUploadedPaths[groupKey] = [];
            for (const file of groups[i].files) {
                const idx = consoleRef.current!.pushLog(`${file.name} ${bar(0)}`);
                const result = await uploadFile(file, idx);
                if (!result.success) {
                    consoleRef.current!.pushLog(`ERROR: Failed to upload ${file.name}`, idx);
                    setIsProcessing(false);
                    return;
                }
                groupUploadedPaths[groupKey].push(result.path);
            }
        }

        const uuid = crypto.randomUUID?.() ?? Math.random().toString(36).slice(2) + Date.now().toString(36);
        const outputFilename = `${uuid}.csv`;

        const params = new URLSearchParams();
        params.append('script', 'spatial-score/003_stats.py');
        for (let i = 0; i < groups.length; i++) {
            const groupKey = `group${i + 1}`;
            params.append(`groupname${i + 1}`, groups[i].name || groupKey);
            for (const path of groupUploadedPaths[groupKey])
                params.append(groupKey, path);
        }
        params.append('output', `data/${outputFilename}`);

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
                if (data.exitCode === 0) setDownloadFilename(outputFilename);
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
                            Spatial Score: Stats
                        </h1>
                        {isProcessing && (
                            <span className="text-xs text-red-400 font-medium uppercase tracking-wider animate-pulse">
                                Do not close this page
                            </span>
                        )}
                    </div>
                    <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
                </header>

                <div className="mb-6 space-y-3">
                    {groups.map((grp, i) => (
                        <div key={i} className="flex items-center gap-3">
                            <input
                                type="text"
                                value={grp.name}
                                onChange={e => updateName(i, e.target.value)}
                                placeholder={`Group ${i + 1} name`}
                                className="px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors w-48"
                            />
                            <input
                                type="file"
                                multiple
                                accept=".csv"
                                onChange={e => updateFiles(i, e.target.files)}
                                className="flex-1 px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                            />
                            {groups.length > MIN_GROUPS && (
                                <button onClick={() => removeGroup(i)}
                                    className="px-2 py-1 bg-slate-800 text-slate-400 text-sm hover:bg-slate-700 transition-colors">
                                    -
                                </button>
                            )}
                        </div>
                    ))}
                    <div className="flex gap-3 items-center">
                        {groups.length < MAX_GROUPS && (
                            <button onClick={addGroup}
                                className="px-3 py-1 bg-slate-800 text-slate-300 text-sm hover:bg-slate-700 transition-colors">
                                +
                            </button>
                        )}
                        <button onClick={handleRun} disabled={isProcessing || !runValid}
                            className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors">
                            {isProcessing ? "Running..." : "Run"}
                        </button>
                        {downloadFilename && !isProcessing && (
                            <a href={`/api/download-file?filename=${downloadFilename}`} download
                                className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors">
                                ↓ Download CSV
                            </a>
                        )}
                    </div>
                </div>

                <Console ref={consoleRef} />
            </div>
        </div>
    );
}
