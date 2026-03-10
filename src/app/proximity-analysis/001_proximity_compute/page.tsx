"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useRef, useState } from "react";

export default function ProximityComputePage() {
    const [headers, setHeaders] = useState<string[][]>([]);
    const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
    const [headerSet, setHeaderSet] = useState<Set<string>>(new Set());
    const [centroidX, setCentroidX] = useState("");
    const [centroidY, setCentroidY] = useState("");
    const [parentArea, setParentArea] = useState("");
    const [parentRegion, setParentRegion] = useState("");
    const [cellType, setCellType] = useState("");
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
        <select value={value} onChange={e => onChange(e.target.value)}>
            <option value="">-- {label} --</option>
            {[...headerSet].map(h => <option key={h} value={h}>{h}</option>)}
        </select>
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
        const outputFilename = `data/${crypto.randomUUID()}.zip`;
        const params = new URLSearchParams({
            script: "proximity-analysis/001_compute.py",
            inputs: uploadedFilenames.join(","),
            centroid_x_col: centroidX,
            centroid_y_col: centroidY,
            parent_area_col: parentArea,
            parent_region_col: parentRegion,
            cell_type_col: cellType,
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
        <div>
            <h1>Proximity Compute</h1>
            <div>
                <input type="file" multiple accept=".csv" onChange={handleFileChange} />
                {/*<pre>{JSON.stringify(headers, null, 2)}</pre>*/}
            </div>
            <div>
                {/* Select centroid x column */}
                {columnSelect(centroidX, setCentroidX, "Select centroid x column")}
            </div>
            <div>
                {/* Select centroid y column */}
                {columnSelect(centroidY, setCentroidY, "Select centroid y column")}
            </div>
            <div>
                {/* Select parent region column */}
                {columnSelect(parentRegion, setParentRegion, "Select parent region column")}
            </div>
            <div>
                {/* Select parent area column */}
                {columnSelect(parentArea, setParentArea, "Select parent area column")}
            </div>
            <div>
                {/* Select cell type column */}
                {columnSelect(cellType, setCellType, "Select cell type column")}
            </div>
            <button onClick={handleRun}>Run</button>
            {downloadFilename && <a href={`/api/download-file?filename=${downloadFilename}`} download>Download</a>}


            <Console ref={consoleRef} />
        </div>
    );
}
