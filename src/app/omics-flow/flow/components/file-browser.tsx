'use client';

import { useEffect } from 'react';
import { useFileBrowser } from '@/hooks/use-file-browser';
import { FileItem } from '../lib/file-browser-session';

export interface FileBrowserProps {
    node_id: string;
    linux_user: string;
    private_key: string;
    host?: string;
    port?: number;
    wsUrl?: string;
    initialPath?: string;
    environment?: string;
}

export default function FileBrowser({
    node_id,
    linux_user,
    private_key,
    host,
    port,
    wsUrl,
    initialPath,
    environment
}: FileBrowserProps) {
    const { connect, toggleDirectory, isConnected, fileItems, error } = useFileBrowser();

    useEffect(() => {
        connect({
            node_id,
            linux_user,
            private_key,
            host,
            port,
            wsUrl,
            initialPath,
            environment
        });
    }, [node_id, linux_user, private_key, host, port, wsUrl, initialPath, environment, connect]);

    const handleToggle = (item: FileItem) => {
        if (item.type === 'directory') {
            toggleDirectory(item);
        }
    };

    const handleDownload = (item: FileItem) => {
        console.log('Download placeholder for:', item.fullPath);
    };

    return (
        <div className="w-full h-full bg-gray-900 text-gray-100 p-4 overflow-auto font-mono text-sm">
            <div className="mb-2 text-xs text-gray-400">
                {isConnected ? '● Connected' : '○ Disconnected'}
                {error && <span className="ml-4 text-red-400">Error: {error}</span>}
            </div>

            <div className="space-y-1">
                {fileItems.map((item, index) => (
                    <div
                        key={`${item.fullPath}-${index}`}
                        className="flex items-center hover:bg-gray-800 px-2 py-1 rounded"
                        style={{ paddingLeft: `${item.depth * 20 + 8}px` }}
                    >
                        {item.type === 'directory' ? (
                            <>
                                <button
                                    onClick={() => handleToggle(item)}
                                    className="mr-2 text-gray-400 hover:text-white w-4 flex-shrink-0"
                                >
                                    {item.isExpanded ? '[-]' : '[+]'}
                                </button>
                                <span className="text-blue-400">📁 {item.name}</span>
                            </>
                        ) : (
                            <>
                                <span className="mr-2 w-4 flex-shrink-0"></span>
                                <span className="text-gray-300 flex-1">📄 {item.name}</span>
                                <button
                                    onClick={() => handleDownload(item)}
                                    className="ml-4 px-2 py-0.5 text-xs bg-gray-700 hover:bg-gray-600 rounded text-gray-300"
                                >
                                    Download
                                </button>
                            </>
                        )}
                    </div>
                ))}
            </div>

            {fileItems.length === 0 && isConnected && (
                <div className="text-gray-500 text-center mt-8">Loading...</div>
            )}
        </div>
    );
}