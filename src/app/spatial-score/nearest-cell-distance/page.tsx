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
    const [columnList, setColumnList] = useState<string[]>([]);

    const [logs, setLogs] = useState<string[]>([]);

    const handleFile = (id: number, e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        setFileItems((prev) => {
            const exists = prev.find((item) => item.id === id);
            if (exists) return prev.map((item) => item.id === id ? { ...item, file } : item);
            return [...prev, { id, file }];
        });
        const reader = new FileReader();
        reader.onload = () => {
            const firstLine = (reader.result as string).split("\n")[0];
            setColumnList(firstLine.split(",").map(c => c.trim()));
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
        <button onClick={handleSubmit}>Submit</button>
        <pre>{JSON.stringify(columnList)}</pre>
        <pre>{logs.join("\n")}</pre>
    </>);
}