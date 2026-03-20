"use client";
import { useState } from "react";

export interface CellData {
    id: string;
    x: number;
    y: number;
    cluster: number;
}

export interface ParsedCsvData {
    data: CellData[];
    dataLines: string[];
    csvHeader: string;
    fileName: string;
    initialMask: Float32Array | null;
}

interface Props {
    onLoad: (parsed: ParsedCsvData) => void;
}

export default function CsvUploadDialog({ onLoad }: Props) {
    const [headers, setHeaders] = useState<string[]>([]);
    const [rawLines, setRawLines] = useState<string[]>([]);
    const [csvHeader, setCsvHeader] = useState("");
    const [xCol, setXCol] = useState("");
    const [yCol, setYCol] = useState("");
    const [clusterCol, setClusterCol] = useState("");
    const [idCol, setIdCol] = useState("");
    const [fileName, setFileName] = useState("");
    const [maskCol, setMaskCol] = useState("");

    const handleFile = async (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        const text = await file.text();
        const lines = text.split(/\r?\n/);
        const header = lines[0];
        const dataLines = lines.slice(1).filter(Boolean);
        setCsvHeader(header);
        setRawLines(dataLines);
        setFileName(file.name);
        setHeaders(header.split(","));
        setXCol(""); setYCol(""); setClusterCol(""); setIdCol(""); setMaskCol("");
    };

    const canLoad = headers.length > 0 && xCol && yCol && clusterCol;

    const handleLoad = () => {
        const xIdx = headers.indexOf(xCol);
        const yIdx = headers.indexOf(yCol);
        const cIdx = headers.indexOf(clusterCol);
        const idIdx   = idCol   ? headers.indexOf(idCol)   : -1;
        const maskIdx = maskCol ? headers.indexOf(maskCol) : -1;
        const data: CellData[] = rawLines.map((line, i) => {
            const cols = line.split(",");
            return {
                id: idIdx >= 0 ? cols[idIdx] : String(i),
                x: parseFloat(cols[xIdx]),
                y: parseFloat(cols[yIdx]),
                cluster: parseInt(cols[cIdx]),
            };
        });
        const initialMask = maskIdx >= 0
            ? new Float32Array(rawLines.map(line => parseFloat(line.split(",")[maskIdx]) || 0))
            : null;
        onLoad({ data, dataLines: rawLines, csvHeader, fileName, initialMask });
    };

    const colSelect = (value: string, onChange: (v: string) => void, label: string) => (
        <div>
            <label style={{ fontSize: 11, color: "#94a3b8", textTransform: "uppercase", letterSpacing: "0.05em", display: "block", marginBottom: 4 }}>{label}</label>
            <select value={value} onChange={e => onChange(e.target.value)}
                style={{ width: "100%", padding: "6px 10px", background: "#0f172a", color: "#e2e8f0", fontSize: 13, border: "none", borderBottom: "1px solid #334155", outline: "none" }}>
                <option value="">— select —</option>
                {headers.map(h => <option key={h} value={h}>{h}</option>)}
            </select>
        </div>
    );

    return (
        <div style={{ position: "fixed", inset: 0, background: "rgba(0,0,0,0.8)", display: "flex", alignItems: "center", justifyContent: "center", zIndex: 1000 }}>
            <div style={{ background: "#1e293b", color: "#e2e8f0", padding: 32, width: 460, fontFamily: "sans-serif", boxShadow: "0 25px 50px rgba(0,0,0,0.5)" }}>
                <h2 style={{ margin: "0 0 4px", fontWeight: 300, fontSize: 20, letterSpacing: "0.05em" }}>Cell Annotator</h2>
                <p style={{ margin: "0 0 24px", fontSize: 13, color: "#94a3b8" }}>Upload a CSV file to begin annotating your spatial data.</p>

                <div style={{ marginBottom: 20 }}>
                    <label style={{ fontSize: 11, color: "#94a3b8", textTransform: "uppercase", letterSpacing: "0.05em", display: "block", marginBottom: 4 }}>CSV File</label>
                    <input type="file" accept=".csv" onChange={handleFile}
                        style={{ width: "100%", padding: "6px 10px", background: "#0f172a", color: "#e2e8f0", fontSize: 13, border: "none", borderBottom: "1px solid #334155", outline: "none", boxSizing: "border-box",
                            // file button styles are set via CSS classes; inline styles don't reach ::file-selector-button
                        }} />
                </div>

                {headers.length > 0 && (
                    <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16, marginBottom: 24 }}>
                        {colSelect(xCol, setXCol, "X Column")}
                        {colSelect(yCol, setYCol, "Y Column")}
                        {colSelect(clusterCol, setClusterCol, "Cluster Column")}
                        {colSelect(idCol, setIdCol, "ID Column (optional)")}
                        {colSelect(maskCol, setMaskCol, "Initial Selection Mask (optional)")}
                    </div>
                )}

                <button onClick={handleLoad} disabled={!canLoad}
                    style={{ padding: "8px 20px", background: canLoad ? "#334155" : "#1e293b", color: canLoad ? "#e2e8f0" : "#475569", fontSize: 13, border: "1px solid #334155", cursor: canLoad ? "pointer" : "not-allowed", transition: "background 0.15s" }}>
                    Load
                </button>
            </div>
        </div>
    );
}
