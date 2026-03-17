"use client";

import { useState, useEffect } from "react";


export default function ComputeRatios() {


    const [headers, setHeaders] = useState<string[][]>([]);
    const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
    const [headerSet, setHeaderSet] = useState<Set<string>>(new Set());

    const [phenotypeColumn, setPhenotypeColumn] = useState("");
    const [imageidColumn, setImageidColumn] = useState("");
    const [phenotypeSet, setPhenotypeSet] = useState<Set<string>>(new Set());

    useEffect(() => {
        if (!phenotypeColumn || !imageidColumn) return;
        const excluded = new Set([phenotypeColumn, imageidColumn]);
        const union = headers.reduce((acc, hdrs) => {
            hdrs.filter(h => !excluded.has(h)).forEach(h => acc.add(h));
            return acc;
        }, new Set<string>());
        setPhenotypeSet(union);
    }, [phenotypeColumn, imageidColumn, headers]);

    const handleFileChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
        const files = Array.from(event.target.files ?? []);
        setSelectedFiles(files);
        const results = await Promise.all(
            files.map(file => file.text().then(text => text.split(/\r?\n/)[0].split(",")))
        );
        setHeaders(results);
        const intersection = results.reduce((acc, cur) => acc.filter(h => cur.includes(h)));
        setHeaderSet(new Set(intersection));
    };

    const columnSelect = (value: string, onChange: (v: string) => void, label: string) => (
        <div>
            <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">{label}</label>
            <select value={value} onChange={e => onChange(e.target.value)}
                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors">
                <option value="">— select —</option>
                {[...headerSet].map(h => <option key={h} value={h}>{h}</option>)}
            </select>
        </div>
    );

    const stateCheck = phenotypeColumn && imageidColumn;

    return (
        <div>
            <div>
                <label className="text-xs text-slate-400 uppercase tracking-wider block mb-1.5">Input Files</label>
                <input type="file" multiple accept=".csv" onChange={handleFileChange}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700" />
            </div>
            <div>
                {columnSelect(phenotypeColumn, setPhenotypeColumn, "Phenotype Column")}
                {columnSelect(imageidColumn, setImageidColumn, "Image ID Column")}
            </div>
            {stateCheck && (
                <div>
                    sdf
                </div>
            )}
        </div>
    )
}