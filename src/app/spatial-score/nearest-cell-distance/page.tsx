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

    const handleFile = (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = () => {
            const firstLine = (reader.result as string).split("\n")[0];
            setColumnList(firstLine.split(",").map(c => c.trim()));
        };
        reader.readAsText(file);
    };

    const handleSubmit = async() => {
        // Step 1: upload the file
    }

    return (<>
        <input type="file" accept=".csv" onChange={handleFile} />
        <pre>{JSON.stringify(columnList)}</pre>
    </>);
}