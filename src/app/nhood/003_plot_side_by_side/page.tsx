import { useRef, useState } from "react";
import { LogDisplay, LogDisplayHandle } from "../components/log-display";


export default function NHoodPlotSideBySide() {
    const [file, setFile] = useState<File | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [outputFilename, setOutputFilename] = useState<string | null>(null);
    const logRef = useRef<LogDisplayHandle>(null);

    return (


                <div className="grid grid-cols-2 gap-6">
                    <LogDisplay ref={logRef} />

                    <div>
                        <div className="flex justify-between items-center mb-3">
                            <h2 className="text-sm font-light text-slate-400 uppercase tracking-wider">
                                Output File
                            </h2>
                            {!isProcessing && outputFilename && (
                                <a
                                    href={`/api/download-file?filename=${outputFilename}`}
                                    download
                                    className="px-3 py-1 bg-emerald-900/50 text-emerald-400 text-xs hover:bg-emerald-900/70 transition-colors"
                                >
                                    ↓ Download CSV
                                </a>
                            )}
                        </div>
                        <div className="bg-slate-900/50 backdrop-blur border-l border-slate-800 h-[calc(100vh-480px)] overflow-auto flex items-center justify-center p-3">
                            {!isProcessing && outputFilename ? (
                                <div className="text-slate-300 text-center">
                                    <div className="text-emerald-400 text-lg mb-2">✓</div>
                                    <div className="text-sm">File ready for download</div>
                                    <div className="text-xs text-slate-500 mt-1">{outputFilename}</div>
                                </div>
                            ) : (
                                <div className="text-slate-600 italic">No file generated yet...</div>
                            )}
                        </div>
                    </div>
                </div>
    );
}