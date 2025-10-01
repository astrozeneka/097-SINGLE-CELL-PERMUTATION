"use client";

import { useState, FormEvent } from "react";

interface FileItem {
  id: number;
  file: File | null;
}

export default function AttractionAvoidanceHeatmap() {
  const [fileItems, setFileItems] = useState<FileItem[]>([
    { id: 1, file: null }
  ]);
  const [nextId, setNextId] = useState(2);
  const [logs, setLogs] = useState<string[]>([]);
  const [isProcessing, setIsProcessing] = useState(false);
  const [zipFilename, setZipFilename] = useState<string | null>(null);

  const addFileItem = () => {
    setFileItems([...fileItems, { id: nextId, file: null }]);
    setNextId(nextId + 1);
  };

  const removeFileItem = (id: number) => {
    if (fileItems.length > 1) {
      setFileItems(fileItems.filter(item => item.id !== id));
    }
  };

  const updateFile = (id: number, file: File | null) => {
    setFileItems(fileItems.map(item =>
      item.id === id ? { ...item, file } : item
    ));
  };

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();

    // Validate all files
    for (const item of fileItems) {
      if (!item.file) {
        alert("Please select all CSV files");
        return;
      }
    }

    setIsProcessing(true);
    setLogs([]);

    try {
      // Upload all files
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

      setLogs((prev) => [...prev, "Starting attraction-avoidance heatmap generation..."]);

      // Run analysis with SSE
      const params = new URLSearchParams({
        filenames: JSON.stringify(uploadedFilenames)
      });
      const eventSource = new EventSource(`/api/run-attraction-avoidance-heatmap?${params.toString()}`);

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
              Attraction-Avoidance Heatmap
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
              CSV Files
            </label>

            {fileItems.map((item, index) => (
              <div key={item.id} className="grid grid-cols-[1fr_auto] gap-3 mb-3">
                <div>
                  {index === 0 && (
                    <label className="block text-xs text-slate-500 mb-1">CSV File</label>
                  )}
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => updateFile(item.id, e.target.files?.[0] || null)}
                    className="w-full px-3 py-1.5 bg-slate-900 text-slate-200 text-sm border-b border-slate-700 focus:border-slate-500 focus:outline-none transition-colors file:mr-3 file:py-1 file:px-2 file:rounded file:border-0 file:text-xs file:bg-slate-800 file:text-slate-300 hover:file:bg-slate-700"
                  />
                </div>

                <div className="self-end">
                  <button
                    type="button"
                    onClick={() => removeFileItem(item.id)}
                    disabled={fileItems.length === 1}
                    className="px-3 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 disabled:bg-slate-900 disabled:text-slate-600 disabled:cursor-not-allowed transition-colors"
                  >
                    Delete
                  </button>
                </div>
              </div>
            ))}

            <button
              type="button"
              onClick={addFileItem}
              className="px-3 py-1.5 bg-slate-800 text-slate-200 text-sm hover:bg-slate-700 transition-colors"
            >
              Add
            </button>
          </div>

          <div className="flex justify-end">
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