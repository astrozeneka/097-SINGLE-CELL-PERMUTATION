import { useEffect, useRef, useState, useCallback } from 'react';
import { FileBrowserSession, FileBrowserConfig, FileItem } from '@/app/omics-flow/flow/lib/file-browser-session';

export function useFileBrowser() {
    const sessionRef = useRef<FileBrowserSession | null>(null);
    const [isConnected, setIsConnected] = useState(false);
    const [fileItems, setFileItems] = useState<FileItem[]>([]);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        sessionRef.current = new FileBrowserSession({
            onConnect: () => {
                setIsConnected(true);
                setError(null);
            },
            onClose: () => {
                setIsConnected(false);
            },
            onError: (message) => {
                setError(message);
            },
            onFileStructureChange: (items) => {
                setFileItems(items);
            }
        });

        return () => {
            sessionRef.current?.close();
        };
    }, []);

    const connect = useCallback((config: FileBrowserConfig) => {
        sessionRef.current?.connect(config);
    }, []);

    const toggleDirectory = useCallback((item: FileItem) => {
        sessionRef.current?.toggleDirectory(item);
    }, []);

    const close = useCallback(() => {
        sessionRef.current?.close();
    }, []);

    return {
        connect,
        toggleDirectory,
        close,
        isConnected,
        fileItems,
        error
    };
}
