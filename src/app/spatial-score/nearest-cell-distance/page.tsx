"use client";

import { useState } from "react";

interface FileItem {
  id: number;
  file: File | null;
}

export default function NearestCellDistance() {
    // Store the input file to be run by the algorithm
    const [fileItems, setFileItems] = useState<FileItem[]>([]);

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
    const [classificationColumn, setClassificationColumn] = useState<string|null>(null);

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
        const uploadedFilenames: string[] = [];

        for (const item of fileItems) {
            setLogs((prev) => [...prev, `Uploading ${item.file!.name}...`]);

            const formData = new FormData();
            formData.append("file", item.file!);
            const uploadResponse = await fetch("/api/upload-file", { method: "POST", body: formData });
            const result = await uploadResponse.json();
            if (!result.success) throw new Error("File upload failed");

            uploadedFilenames.push(result.filename);
            setLogs((prev) => [...prev, `Uploaded ${item.file!.name} as ${result.filename}`]);
        }
    }

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
        <select disabled={!columnList.size} value={classificationColumn ?? ""} onChange={(e) => setClassificationColumn(e.target.value)}>
            <option value="">Classification</option>
            {[...columnList].map(c => <option key={c} value={c}>{c}</option>)}
        </select>
        <hr/>
        <button onClick={handleSubmit}>Submit</button>
        <pre>{logs.join("\n")}</pre>
    </>);
}