"use client";

import { useEffect, useRef, useState } from "react";

interface FileItem {
  id: number;
  file: File | null;
}

export default function NearestCellDistance() {
    // Store the input file to be run by the algorithm
    const [fileItems, setFileItems] = useState<FileItem[]>([]);
    const [isProcessing, setIsProcessing] = useState(false);
    const eventSourceRef = useRef<EventSource | null>(null);

    // The list of columns
    const [columnList, setColumnList] = useState<Set<string>>(new Set());

    // Logs
    const [logs, setLogs] = useState<string[]>([]);

    // Error message for column mismatch
    const [error, setError] = useState<string | null>(null);

    // Columns to be selected
    const [objectIdColumn, setObjectIdColumn] = useState<string|null>(null);
    const [centroidXColumn, setCentroidXColumn] = useState<string|null>(null);
    const [centroidYColumn, setCentroidYColumn] = useState<string|null>(null);
    const [phenotypeColumn, setPhenotypeColumn] = useState<string|null>(null);
    const [imageIdColumn, setImageIdColumn] = useState<string|null>(null);

    const getDataColumns = (text: string): Set<string> =>
        new Set(text.split("\n")[0].split(",").map(c => c.trim()));

    const handleFile = (id: number, e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        setError(null);
        const reader = new FileReader();
        reader.onload = () => {
            const cols = getDataColumns(reader.result as string);
            setFileItems((prev) => {
                const exists = prev.find((item) => item.id === id);
                if (exists) return prev.map((item) => item.id === id ? { ...item, file } : item);
                return [...prev, { id, file }];
            });
            setColumnList((prev) => prev.size === 0 ? cols : new Set([...prev].filter(c => cols.has(c))));
        };
        reader.readAsText(file);
    };

    const handleSubmit = async () => {
        setIsProcessing(true);
        setLogs([]);

        try {
            const uploadedFilenames: string[] = [];

            for (const item of fileItems) {
                setLogs((prev) => [...prev, `Uploading ${item.file!.name}...`]);

                const formData = new FormData();
                formData.append("file", item.file!);
                const uploadResponse = await fetch("/api/upload-file", { method: "POST", body: formData });
                const result = await uploadResponse.json();
                if (!result.success) throw new Error("File upload failed");

                uploadedFilenames.push(result.filename);
            }

            setLogs((prev) => [...prev, "Starting script..."]);
            const params = new URLSearchParams({
                script: "spatial_score/001_compute_distances.py",
                ...Object.fromEntries(uploadedFilenames.map((f, i) => [`input${i + 1}`, `data/${f}`])),
                "object-id": objectIdColumn!,
                "centroid-x": centroidXColumn!,
                "centroid-y": centroidYColumn!,
                phenotype: phenotypeColumn!,
                imageid: imageIdColumn!,
            });

            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
            eventSourceRef.current = eventSource;

            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);
                if (data.type === "stdout") {
                    setLogs((prev) => [...prev, data.content]);
                } else if (data.type === "stderr") {
                    setLogs((prev) => [...prev, `ERROR: ${data.content}`]);
                } else if (data.type === "complete") {
                    setLogs((prev) => [...prev, `Process completed with exit code: ${data.exitCode}`]);
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    setLogs((prev) => [...prev, `ERROR: ${data.message}`]);
                    eventSource.close();
                    eventSourceRef.current = null;
                    setIsProcessing(false);
                }
            };
            eventSource.onerror = () => {
                setLogs((prev) => [...prev, "Connection error"]);
                eventSource.close();
                eventSourceRef.current = null;
                setIsProcessing(false);
            };
        } catch (error: any) {
            setLogs((prev) => [...prev, `Error: ${error.message}`]);
            setIsProcessing(false);
        }
    };

    useEffect(() => {
        return () => { eventSourceRef.current?.close(); };
    }, []);

    return (<>
        <input type="file" accept=".csv" onChange={(e) => handleFile(1, e)} />
        <input type="file" accept=".csv" onChange={(e) => handleFile(2, e)} />
        <input type="file" accept=".csv" onChange={(e) => handleFile(3, e)} />
        <input type="file" accept=".csv" onChange={(e) => handleFile(4, e)} />
        {error && <p style={{ color: "red" }}>{error}</p>}
        <hr/>
        <select disabled={!columnList.size} value={objectIdColumn ?? ""} onChange={(e) => setObjectIdColumn(e.target.value)}>
            <option value="">Object ID</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <select disabled={!columnList.size} value={centroidXColumn ?? ""} onChange={(e) => setCentroidXColumn(e.target.value)}>
            <option value="">Centroid X</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <select disabled={!columnList.size} value={centroidYColumn ?? ""} onChange={(e) => setCentroidYColumn(e.target.value)}>
            <option value="">Centroid Y</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <select disabled={!columnList.size} value={phenotypeColumn ?? ""} onChange={(e) => setPhenotypeColumn(e.target.value)}>
            <option value="">Phenotype</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <select disabled={!columnList.size} value={imageIdColumn ?? ""} onChange={(e) => setImageIdColumn(e.target.value)}>
            <option value="">Image ID</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <hr/>
        <button onClick={handleSubmit} disabled={isProcessing || !fileItems.length}>
            {isProcessing ? "Running..." : "Submit"}
        </button>
        <pre>{logs.join("\n")}</pre>
    </>);
}