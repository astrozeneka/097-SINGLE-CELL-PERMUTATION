"use client";

import { useState, FormEvent } from "react";

export default function Home() {
  const [file, setFile] = useState<File | null>(null);
  const [nPermutations, setNPermutations] = useState<number>(10);
  const [logs, setLogs] = useState<string[]>([]);
  const [isProcessing, setIsProcessing] = useState(false);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    if (!file) {
      alert("Please select a file");
      return;
    }

    setIsProcessing(true);
    setLogs([]);

    try {
      // Upload file
      const formData = new FormData();
      formData.append("file", file);

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

      // Run permutation test with SSE
      setLogs((prev) => [...prev, "Starting permutation test..."]);
      const eventSource = new EventSource("/api/run-permutation");

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
    <div className="min-h-screen bg-gray-50 p-8">
      <div className="max-w-4xl mx-auto">
        <h1 className="text-3xl font-bold text-gray-800 mb-8">
          Spatial Distance Permutation Test
        </h1>

        <form onSubmit={handleSubmit} className="bg-white rounded-lg shadow-md p-6 mb-8">
          <div className="mb-4">
            <label className="block text-gray-700 font-medium mb-2">
              File:
            </label>
            <input
              type="file"
              accept=".csv"
              onChange={(e) => setFile(e.target.files?.[0] || null)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
          </div>

          <div className="mb-6">
            <label className="block text-gray-700 font-medium mb-2">
              N Permutations:
            </label>
            <input
              type="number"
              value={nPermutations}
              onChange={(e) => setNPermutations(parseInt(e.target.value))}
              min="1"
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
          </div>

          <button
            type="submit"
            disabled={isProcessing}
            className="w-full bg-blue-600 text-white font-medium py-2 px-4 rounded-md hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed transition-colors"
          >
            {isProcessing ? "Processing..." : "Submit"}
          </button>
        </form>

        <div className="bg-white rounded-lg shadow-md p-6">
          <h2 className="text-xl font-bold text-gray-800 mb-4">Output Logs</h2>
          <div className="bg-gray-900 text-green-400 font-mono text-sm p-4 rounded-md h-96 overflow-y-auto">
            {logs.length === 0 ? (
              <div className="text-gray-500">No logs yet...</div>
            ) : (
              logs.map((log, index) => (
                <div key={index} className="mb-1">
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
