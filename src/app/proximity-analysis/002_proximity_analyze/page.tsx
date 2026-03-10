"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useEffect, useRef, useState } from "react";


export default function ProximityComputePage() {

    const [group1Name, setGroup1Name] = useState("Group 1");
    const [group2Name, setGroup2Name] = useState("Group 2");
    const [group1Files, setGroup1Files] = useState<File[]>([]);
    const [group2Files, setGroup2Files] = useState<File[]>([]);
    const consoleRef = useRef<ConsoleHandle>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const eventSourceRef = useRef<EventSource | null>(null);
    const [downloadCSVFilename, setDownloadCSVFilename] = useState<string | null>(null);

    const handleFilesChange = (event: React.ChangeEvent<HTMLInputElement>, mutator: React.Dispatch<React.SetStateAction<File[]>>) => {
        const files = Array.from(event.target.files ?? []);
        mutator(files);
        // No need to compute here
    };

    const handleGroup1FilesChange = (e: React.ChangeEvent<HTMLInputElement>) => handleFilesChange(e, setGroup1Files);
    const handleGroup2FilesChange = (e: React.ChangeEvent<HTMLInputElement>) => handleFilesChange(e, setGroup2Files);


    // When the downloadCSVFilename is set
    useEffect(() => {
        if (!downloadCSVFilename) return;
        
        // Process the filename to create a download link
        console.log("here")
    }, [downloadCSVFilename]);
    
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
        
            const uploadedFilePaths: string[] = [];
            for (const files of [group1Files, group2Files]) {
                for (const file of files) {
                    const idx = consoleRef.current!.pushLog(`${file.name} ${bar(0)}`);
                    const uploadResult = await uploadFile(file, idx);
                    if (!uploadResult.success) {
                        consoleRef.current!.pushLog(`ERROR: Failed to upload ${file.name}`, idx);
                        setIsProcessing(false);
                        return;
                    }
                    uploadedFilePaths.push(uploadResult.path);
                }
            }

            // Pass relevant parameters
            const outputFilename = `data/${crypto.randomUUID()}.csv`;
            const params = new URLSearchParams({
                script: "proximity-analysis/002_analyze.py",
                group1_inputs: group1Files.map(f => `data/${f.name}`).join(","),
                group2_inputs: group2Files.map(f => `data/${f.name}`).join(","),
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
                <input type="file" multiple accept=".csv" onChange={handleGroup1FilesChange} />
                <input type="file" multiple accept=".csv" onChange={handleGroup2FilesChange} />
                {/*<pre>{JSON.stringify(headers, null, 2)}</pre>*/}
            </div>

            <button onClick={handleRun}>Run</button>
            <table>
                <thead>
                    <tr>
                        <th>Cell Type A</th>
                        <th>Cell Type B</th>
                        <th>Interaction/mm2 in {group1Name}</th>
                        <th>Interaction/mm2 in {group2Name}</th>
                        <th>p-value</th>
                    </tr>
                </thead>
                <tbody>

                </tbody>
            </table>
            
            <Console ref={consoleRef} />
        </div>
    );
}