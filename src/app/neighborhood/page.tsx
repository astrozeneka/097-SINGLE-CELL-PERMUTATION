
"use client";
import { useState } from "react";


export default function NeighborhoodPage() {


    const [file, setFile] = useState<File | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [logs, setLogs] = useState<string[]>([]);

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();

        try {
            // Upload file
            const formData = new FormData();
            formData.append("file", file as any);

            setLogs((prev) => [...prev, "Uploading file..."]);
            const uploadResponse = await fetch("/api/upload-file", {
                method: "POST",
                body: formData,
            });

            const uploadResult = await uploadResponse.json();
            if (!uploadResult.success) {
                throw new Error("File upload failed");
            }

            setLogs((prev) => [...prev, `File uploaded: ${uploadResult.filename}`]);

            // Run neighborhood analysis with SSE
            setLogs((prev) => [...prev, "Starting neighborhood analysis..."]);
            let output_name = uploadResult.filename.replace('.csv', '_cnhood.csv');
            const params = new URLSearchParams({
                input: uploadResult.filename,
                output: output_name,
                script: "nhood/001_neighborhood.py"
                // TODO later, parameter, radius/knn
            });
            const eventSource = new EventSource(`/api/run-python?${params.toString()}`);
            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);
                console.log(data)
            };

        } catch (error: any) {
            // TODO
            setLogs((prev) => [...prev, `Error: ${error.message}`]);
            console.error("Error:", error);
        }

    }

    return (
        <form onSubmit={handleSubmit} className="mb-6 space-y-4">

            <div>
                <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                    Input File
                </label>
                <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => setFile(e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                />
            </div>


            <button
                type="submit"
                disabled={isProcessing}
                className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
            >
                {isProcessing ? "Running..." : "Generate"}
            </button>

        </form>
    )
}