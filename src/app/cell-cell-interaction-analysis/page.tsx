"use client";

import { useState } from "react";

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

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100 p-6">
      <div className="max-w-5xl mx-auto">
        <header className="mb-6">
          <h1 className="text-2xl font-light text-slate-100 tracking-wide">
            Cell-Cell Interaction Analysis
          </h1>
          <div className="h-px bg-gradient-to-r from-slate-700 via-slate-600 to-transparent mt-2"></div>
        </header>

        <div className="mb-6">
          <label className="block text-xs text-slate-400 mb-3 uppercase tracking-wider">
            File Pairs
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
      </div>
    </div>
  );
}
