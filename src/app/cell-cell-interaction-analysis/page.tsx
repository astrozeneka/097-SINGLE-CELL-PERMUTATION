"use client";

import { useState, FormEvent } from "react";

interface FileTuple {
  id: number;
  file1: File | null;
  file2: File | null;
}

export default function CellCellInteractionAnalysis() {
  const [fileTuples, setFileTuples] = useState<FileTuple[]>([
    { id: 1, file1: null, file2: null }
  ]);
  const [nextId, setNextId] = useState(2);
  const [logs, setLogs] = useState<string[]>([]);
  const [isProcessing, setIsProcessing] = useState(false);
  const [outputDir, setOutputDir] = useState<string>("attraction_repulsion_results");
  const [zipFilename, setZipFilename] = useState<string | null>(null);

  const addFileTuple = () => {
    setFileTuples([...fileTuples, { id: nextId, file1: null, file2: null }]);
    setNextId(nextId + 1);
  };

  const removeFileTuple = (id: number) => {
    if (fileTuples.length > 1) {
      setFileTuples(fileTuples.filter(tuple => tuple.id !== id));
    }
  };

  const updateFile1 = (id: number, file: File | null) => {
    setFileTuples(fileTuples.map(tuple =>
      tuple.id === id ? { ...tuple, file1: file } : tuple
    ));
  };

  const updateFile2 = (id: number, file: File | null) => {
    setFileTuples(fileTuples.map(tuple =>
      tuple.id === id ? { ...tuple, file2: file } : tuple
    ));
  };

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();

    // Validate all file tuples
    for (const tuple of fileTuples) {
      if (!tuple.file1 || !tuple.file2) {
        alert("Please select both files for each pair");
        return;
      }
    }

    setIsProcessing(true);
    setLogs([]);

    try {
      // Upload all files and build file tuple data
      const uploadedFileTuples: Array<[string, string, string]> = [];

      for (const tuple of fileTuples) {
        setLogs((prev) => [...prev, `Uploading files for pair ${tuple.id}...`]);

        // Upload file1 (observed)
        const formData1 = new FormData();
        formData1.append("file", tuple.file1!);
        const upload1 = await fetch("/api/upload-file", { method: "POST", body: formData1 });
        const result1 = await upload1.json();
        if (!result1.success) throw new Error("File upload failed");

        // Upload file2 (permuted)
        const formData2 = new FormData();
        formData2.append("file", tuple.file2!);
        const upload2 = await fetch("/api/upload-file", { method: "POST", body: formData2 });
        const result2 = await upload2.json();
        if (!result2.success) throw new Error("File upload failed");

        // Generate slug from file1 name (remove extension)
        const slug = tuple.file1!.name.replace(/\.[^/.]+$/, "");
        uploadedFileTuples.push([slug, result1.filename, result2.filename]);
      }

      setLogs((prev) => [...prev, "Starting cell-cell interaction analysis..."]);

      // Run analysis with SSE
      const params = new URLSearchParams({
        file_tuples: JSON.stringify(uploadedFileTuples),
        output_dir: outputDir
      });
      const eventSource = new EventSource(`/api/run-cell-cell-interaction?${params.toString()}`);

      eventSource.onmessage = (event) => {
        const data = JSON.parse(event.data);

        if (data.type === "acknowledgment") {
          setLogs((prev) => [...prev, data.message]);
        } else if (data.type === "stdout") {
          setLogs((prev) => [...prev, data.content]);
        } else if (data.type === "stderr") {
          setLogs((prev) => [...prev, `ERROR: ${data.content}`]);
        } else if (data.type === "complete") {
          setLogs((prev) => [...prev, `Process completed with exit code: ${data.exitCode}`]);
          if (data.zipFilename) {
            setZipFilename(data.zipFilename);
          }
          eventSource.close();
          setIsProcessing(false);
        } else if (data.type === "error") {
          setLogs((prev) => [...prev, `ERROR: ${data.message}`]);
          eventSource.close();
          setIsProcessing(false);
        }
      };

      eventSource.onerror = () => {
        setLogs((prev) => [...prev, "Connection error"]);
        eventSource.close();
        setIsProcessing(false);
      };
    } catch (error: any) {
      setLogs((prev) => [...prev, `Error: ${error.message}`]);
      setIsProcessing(false);
    }
  };

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
      <div className="max-w-5xl mx-auto">
        <header className="mb-6">
          <div className="flex items-center gap-4">
            <h1 className="text-2xl font-light text-slate-100 tracking-wide">
              Cell-Cell Interaction Analysis
            </h1>
            {isProcessing && (
              <span className="text-xs text-red-400 font-medium uppercase tracking-wider animate-pulse">
                Do not close this page
              </span>
            )}
          </div>
          <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
        </header>

        <form onSubmit={handleSubmit} className="mb-6">
          <div className="mb-6">
            <label className="block text-xs text-slate-400 mb-3 uppercase tracking-wider">
              File Pairs (Observed / Permuted)
            </label>

            {fileTuples.map((tuple) => (
              <div key={tuple.id} className="grid grid-cols-[1fr_1fr_auto] gap-3 mb-3 items-end">
                <div>
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => updateFile1(tuple.id, e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                  />
                </div>

                <div>
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => updateFile2(tuple.id, e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                  />
                </div>

                <button
                  type="button"
                  onClick={() => removeFileTuple(tuple.id)}
                  disabled={fileTuples.length === 1}
                  className="px-3 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
                >
                  Delete
                </button>
              </div>
            ))}

            <button
              type="button"
              onClick={addFileTuple}
              className="px-3 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 transition-colors"
            >
              Add
            </button>
          </div>

          <div className="grid grid-cols-[1fr_auto] gap-3 items-end">
            <div>
              <label className="block text-xs text-slate-400 mb-1.5 uppercase tracking-wider">
                Output Directory
              </label>
              <input
                type="text"
                value={outputDir}
                onChange={(e) => setOutputDir(e.target.value)}
                className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
              />
            </div>

            <button
              type="submit"
              disabled={isProcessing}
              className="px-5 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
            >
              {isProcessing ? "Running..." : "Run"}
            </button>
          </div>
        </form>

        <div>
          <div className="flex justify-between items-center mb-3">
            <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
              Console Output
            </h2>
            {!isProcessing && zipFilename && (
              <a
                href={`/api/download-file?filename=${zipFilename}`}
                download
                className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors"
              >
                ↓ Download Results
              </a>
            )}
          </div>
          <div className="bg-slate-900/50 backdrop-blur font-mono text-xs p-3 h-[calc(100vh-400px)] overflow-y-auto flex flex-col-reverse border-l border-slate-800">
            {logs.length === 0 ? (
              <div className="text-slate-600 italic">Waiting for input...</div>
            ) : (
              [...logs].reverse().map((log, index) => (
                <div key={index} className="text-slate-300 leading-relaxed py-0.5">
                  {log}
                </div>
              ))
            )}
          </div>
        </div>
      </div>
    </div>
  );
}
