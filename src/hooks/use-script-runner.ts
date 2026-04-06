import { useEffect, useRef, useCallback, useState } from 'react';
import { OMICS_FLOW_WS_URL } from '@/config/api';

export interface ScriptRunnerConfig {
    node_id: string;
    linux_user: string;
    private_key: string;
    environment: string;
    script?: string;
    args?: string[];
    host?: string;
    port?: number;
    wsUrl?: string;
}

export interface ScriptRunnerCallbacks {
    onStdout?: (data: string) => void;
    onStderr?: (data: string) => void;
    onStatus?: (message: string) => void;
    onError?: (message: string) => void;
    onDone?: (exitCode: number) => void;
    onConnect?: () => void;
    onClose?: () => void;
}

export function useScriptRunner(callbacks?: ScriptRunnerCallbacks) {
    const ws = useRef<WebSocket | null>(null);
    const [isConnected, setIsConnected] = useState(false);
    const [isRunning, setIsRunning] = useState(false);

    const run = useCallback((config: ScriptRunnerConfig) => {
        if (ws.current?.readyState === WebSocket.OPEN) {
            ws.current.close();
        }

        const wsUrl = config.wsUrl || `${OMICS_FLOW_WS_URL}/script-runner`;
        console.log("wsUrl:", wsUrl);
        ws.current = new WebSocket(wsUrl);

        ws.current.onopen = () => {
            setIsConnected(true);
            callbacks?.onConnect?.();

            const message = {
                type: 'run',
                node_id: config.node_id,
                linux_user: config.linux_user,
                private_key: config.private_key,
                environment: config.environment,
                script: config.script,
                args: config.args || [],
                host: config.host || '127.0.0.1',
                port: config.port || 22,
            };

            ws.current!.send(JSON.stringify(message));
            setIsRunning(true);
        };

        ws.current.onmessage = (event) => {
            const message = JSON.parse(event.data);

            switch (message.type) {
                case 'stdout':
                    callbacks?.onStdout?.(message.data);
                    break;
                case 'stderr':
                    callbacks?.onStderr?.(message.data);
                    break;
                case 'status':
                    callbacks?.onStatus?.(message.data);
                    break;
                case 'error':
                    callbacks?.onError?.(message.data);
                    setIsRunning(false);
                    break;
                case 'done':
                    callbacks?.onDone?.(message.data.exitCode);
                    setIsRunning(false);
                    break;
            }
        };

        ws.current.onerror = (error) => {
            console.error('WebSocket error:', error);
            callbacks?.onError?.('WebSocket connection error');
            setIsRunning(false);
        };

        ws.current.onclose = () => {
            setIsConnected(false);
            setIsRunning(false);
            callbacks?.onClose?.();
        };
    }, [callbacks]);

    const sendInput = useCallback((data: string) => {
        if (ws.current?.readyState === WebSocket.OPEN) {
            ws.current.send(JSON.stringify({ type: 'input', data }));
        }
    }, []);

    const close = useCallback(() => {
        if (ws.current) {
            ws.current.close();
            ws.current = null;
        }
    }, []);

    useEffect(() => {
        return () => {
            if (ws.current) {
                ws.current.close();
            }
        };
    }, []);

    return {
        run,
        sendInput,
        close,
        isConnected,
        isRunning,
    };
}
