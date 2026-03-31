"use client";
import { useState } from "react";

export interface CellDataV2_1 {
    id: string;
    x: number;
    y: number;
    umap_1: number;
    umap_2: number;
    cluster: string;
    clusterIdx: number;
    SampleId: string;
}

interface Props {
    onLoad: (data: CellDataV2_1[]) => void;
}

export default function CsvUploadDialogV2_1({ onLoad }: Props) {
    const [headers, setHeaders] = useState<string[]>([]);
    const [rawLines, setRawLines] = useState<string[]>([]);
    const [xCol, setXCol] = useState("");
    const [yCol, setYCol] = useState("");
    const [umap1Col, setUmap1Col] = useState("");
    const [umap2Col, setUmap2Col] = useState("");
    const [clusterCol, setClusterCol] = useState("");
    const [sampleIdCol, setSampleIdCol] = useState("");
    const [idCol, setIdCol] = useState("");

    const handleFile = async (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        const text = await file.text();
        const lines = text.split(/\r?\n/);
        const dataLines = lines.slice(1).filter(Boolean);
        setRawLines(dataLines);
        setHeaders(lines[0].split(","));
        setXCol(""); setYCol(""); setUmap1Col(""); setUmap2Col(""); setClusterCol(""); setSampleIdCol(""); setIdCol("");
    };

    const canLoad = headers.length > 0 && xCol && yCol && umap1Col && umap2Col && clusterCol && sampleIdCol;

    const handleLoad = () => {
        const xIdx = headers.indexOf(xCol);
        const yIdx = headers.indexOf(yCol);
        const u1Idx = headers.indexOf(umap1Col);
        const u2Idx = headers.indexOf(umap2Col);
        const cIdx = headers.indexOf(clusterCol);
        const sIdx = headers.indexOf(sampleIdCol);
        const idIdx = idCol ? headers.indexOf(idCol) : -1;

        const clusterOrder: string[] = [];
        const clusterMap = new Map<string, number>();
        const parsed = rawLines.map((line, i) => {
            const cols = line.split(",");
            const cluster = cols[cIdx];
            if (!clusterMap.has(cluster)) { clusterMap.set(cluster, clusterOrder.length); clusterOrder.push(cluster); }
            return {
                id: idIdx >= 0 ? cols[idIdx] : String(i),
                x: +cols[xIdx],
                y: +cols[yIdx],
                umap_1: +cols[u1Idx],
                umap_2: +cols[u2Idx],
                cluster,
                clusterIdx: clusterMap.get(cluster)!,
                SampleId: cols[sIdx],
            };
        });
        onLoad(parsed);
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
                        style={{ width: "100%", padding: "6px 10px", background: "#0f172a", color: "#e2e8f0", fontSize: 13, border: "none", borderBottom: "1px solid #334155", outline: "none", boxSizing: "border-box" }} />
                </div>

                {headers.length > 0 && (
                    <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16, marginBottom: 24 }}>
                        {colSelect(xCol, setXCol, "X Column")}
                        {colSelect(yCol, setYCol, "Y Column")}
                        {colSelect(umap1Col, setUmap1Col, "UMAP 1")}
                        {colSelect(umap2Col, setUmap2Col, "UMAP 2")}
                        {colSelect(clusterCol, setClusterCol, "Cluster Column")}
                        {colSelect(sampleIdCol, setSampleIdCol, "Sample ID Column")}
                        {colSelect(idCol, setIdCol, "ID Column (optional)")}
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
