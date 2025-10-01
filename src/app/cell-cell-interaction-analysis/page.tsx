"use client";

import { useState, FormEvent } from "react";

interface FileTuple {
  id: number;
  sampleId: string;
  observedFile: File | null;
  permutedFile: File | null;
}

export default function CellCellInteractionAnalysis() {
  const [fileTuples, setFileTuples] = useState<FileTuple[]>([
    { id: 1, sampleId: "", observedFile: null, permutedFile: null }
  ]);
  const [nextId, setNextId] = useState(2);
  const [logs, setLogs] = useState<string[]>([]);
  const [isProcessing, setIsProcessing] = useState(false);
  const [outputDir, setOutputDir] = useState<string>("attraction_repulsion_results");
  const [zipFilename, setZipFilename] = useState<string | null>(null);

  const addFileTuple = () => {
    setFileTuples([...fileTuples, { id: nextId, sampleId: "", observedFile: null, permutedFile: null }]);
    setNextId(nextId + 1);
  };

  const removeFileTuple = (id: number) => {
    if (fileTuples.length > 1) {
      setFileTuples(fileTuples.filter(tuple => tuple.id !== id));
    }
  };

  const updateSampleId = (id: number, sampleId: string) => {
    setFileTuples(fileTuples.map(tuple =>
      tuple.id === id ? { ...tuple, sampleId } : tuple
    ));
  };

  const updateObservedFile = (id: number, file: File | null) => {
    setFileTuples(fileTuples.map(tuple =>
      tuple.id === id ? { ...tuple, observedFile: file } : tuple
    ));
  };

  const updatePermutedFile = (id: number, file: File | null) => {
    setFileTuples(fileTuples.map(tuple =>
      tuple.id === id ? { ...tuple, permutedFile: file } : tuple
    ));
  };

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();

    // Validate all file tuples
    for (const tuple of fileTuples) {
      if (!tuple.sampleId.trim()) {
        alert("Please provide a sample ID for each pair");
        return;
      }
      if (!tuple.observedFile || !tuple.permutedFile) {
        alert("Please select both observed and permuted files for each pair");
        return;
      }
    }

    setIsProcessing(true);
    setLogs([]);

    try {
      // Upload all files and build file tuple data
      const uploadedFileTuples: Array<[string, string, string]> = [];

      for (const tuple of fileTuples) {
        setLogs((prev) => [...prev, `Uploading files for sample ${tuple.sampleId}...`]);

        // Upload observed file
        const formDataObserved = new FormData();
        formDataObserved.append("file", tuple.observedFile!);
        const uploadObserved = await fetch("/api/upload-file", { method: "POST", body: formDataObserved });
        const resultObserved = await uploadObserved.json();
        if (!resultObserved.success) throw new Error("Observed file upload failed");

        // Upload permuted file
        const formDataPermuted = new FormData();
        formDataPermuted.append("file", tuple.permutedFile!);
        const uploadPermuted = await fetch("/api/upload-file", { method: "POST", body: formDataPermuted });
        const resultPermuted = await uploadPermuted.json();
        if (!resultPermuted.success) throw new Error("Permuted file upload failed");

        // Format: [sampleId, observedFile, permutedFile]
        uploadedFileTuples.push([tuple.sampleId, resultObserved.filename, resultPermuted.filename]);
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
              Sample Data (Sample ID, Observed File, Permuted File)
            </label>

            {fileTuples.map((tuple, index) => (
              <div key={tuple.id} className="grid grid-cols-[150px_1fr_1fr_auto] gap-3 mb-3">
                <div>
                  {index === 0 && (
                    <label className="block text-xs text-slate-500 mb-1">Sample ID</label>
                  )}
                  <input
                    type="text"
                    value={tuple.sampleId}
                    onChange={(e) => updateSampleId(tuple.id, e.target.value)}
                    placeholder="Sample ID"
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors"
                  />
                </div>

                <div>
                  {index === 0 && (
                    <label className="block text-xs text-slate-500 mb-1">Observed File</label>
                  )}
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => updateObservedFile(tuple.id, e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                  />
                </div>

                <div>
                  {index === 0 && (
                    <label className="block text-xs text-slate-500 mb-1">Permuted File</label>
                  )}
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => updatePermutedFile(tuple.id, e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                  />
                </div>

                <div className="self-end">
                  <button
                    type="button"
                    onClick={() => removeFileTuple(tuple.id)}
                    disabled={fileTuples.length === 1}
                    className="px-3 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
                  >
                    Delete
                  </button>
                </div>
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
