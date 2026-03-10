"use client";

import Console, { ConsoleHandle } from "@/app/components/console";
import { useRef, useState } from "react";

export default function ProximityComputePage() {
    const [headers, setHeaders] = useState<string[][]>([]);
    const [headerSet, setHeaderSet] = useState<Set<string>>(new Set());
    const [centroidX, setCentroidX] = useState("");
    const [centroidY, setCentroidY] = useState("");
    const [parentArea, setParentArea] = useState("");
    const [cellType, setCellType] = useState("");
    const consoleRef = useRef<ConsoleHandle>(null);
    const eventSourceRef = useRef<EventSource | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);

    const handleFileChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
        const files = Array.from(event.target.files ?? []);
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

    const handleRun = () => {
        setIsProcessing(true);
        consoleRef.current?.clearLogs();
        // TODO: pass relevant parameters
        const eventSource = new EventSource("/api/run-python?script=proximity-analysis/compute.py");
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
                {/* Select parent area column */}
                {columnSelect(parentArea, setParentArea, "Select parent area column")}
            </div>
            <div>
                {/* Select cell type column */}
                {columnSelect(cellType, setCellType, "Select cell type column")}
            </div>
            <button onClick={handleRun}>Run</button>


            <Console ref={consoleRef} />
        </div>
    );
}
