"use client";

import { FormEvent, useRef, useState } from "react";
import { LogDisplay, LogDisplayHandle } from "../components/log-display";


export default function NHoodClusterPage() {
    const [file, setFile] = useState<File | null>(null);
    const [isProcessing, setIsProcessing] = useState(false);
    const [outputFilename, setOutputFilename] = useState<string | null>(null);
    const logRef = useRef<LogDisplayHandle>(null);

    // The parameters for the submission
    const [nMotifs, setNMotifs] = useState(20);
    const [knnSeed, setKnnSeed] = useState(42);

    // 
    const handleSubmit = async (e: FormEvent) => {
        e.preventDefault();
        if (!file) {
            alert("Please select a file");
            return;
        }
        setIsProcessing(true);
        logRef.current?.clearLogs();
        setOutputFilename(null);

        try {
            // Upload file
            const formData = new FormData();
            formData.append("file", file);

            logRef.current?.pushLog("Uploading file...");
            const uploadResponse = await fetch("/api/upload-file", {
                method: "POST",
                body: formData,
            });

            const uploadResult = await uploadResponse.json();
            if (!uploadResult.success) {
                setIsProcessing(false);
                logRef.current?.pushLog("File upload failed");
                throw new Error("File upload failed");
            }

            logRef.current?.pushLog("File uploaded successfully");

            // Run the clustering process
            // Begin by preparing the request
            logRef.current?.pushLog("Starting script...");
            const outputName = uploadResult.filename.replace('.csv', '_withmotifs.csv');
            const params = new URLSearchParams({
                input: `data/${uploadResult.filename}`,
                output: `data/${outputName}`,
                n_motifs: nMotifs.toString(),
                knn_seed: knnSeed.toString(),
            });
            const eventSource = new EventSource(`/api/run-nhood-cluster?${params.toString()}`);
            eventSource.onmessage = (event) => {
                const data = JSON.parse(event.data);

                if (data.type === "acknowledgment") {
                    logRef.current?.pushLog(data.message);
                } else if (data.type === "stdout") {
                    logRef.current?.pushLog(data.content);
                } else if (data.type === "stderr") {
                    logRef.current?.pushLog(`ERROR: ${data.content}`);
                } else if (data.type === "complete") {
                    logRef.current?.pushLog(`Process completed with exit code: ${data.exitCode}`);
                    if (data.exitCode === 0 && data.outputFilename) {
                        setOutputFilename(data.outputFilename.replace('data/', ''));
                    }
                    eventSource.close();
                    setIsProcessing(false);
                } else if (data.type === "error") {
                    logRef.current?.pushLog(`ERROR: ${data.message}`);
                    eventSource.close();
                    setIsProcessing(false);
                }
            };
            eventSource.onerror = () => {
                logRef.current?.pushLog("Connection error");
                eventSource.close();
                setIsProcessing(false);
            };

        } catch (error: any) {
            logRef.current?.pushLog(`Error: ${error.message}`);
            setIsProcessing(false);
        } finally {
            setIsProcessing(false);
        }

    }

    return (
        <LogDisplay ref={logRef} />
    )
}