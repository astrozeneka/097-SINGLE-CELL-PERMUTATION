import { OMICS_FLOW_WS_URL } from "@/config/api";

export interface FileBrowserConfig {
    node_id: string;
    linux_user: string;
    private_key: string;
    host?: string;
    port?: number;
    wsUrl?: string;
    initialPath?: string;
    environment?: string;
}

export interface FileItem {
    name: string;
    fullPath: string;
    type: 'file' | 'directory';
    isExpanded: boolean;
    depth: number;
}

export interface FileBrowserCallbacks {
    onConnect?: () => void;
    onClose?: () => void;
    onError?: (message: string) => void;
    onFileStructureChange?: (items: FileItem[]) => void;
}

export class FileBrowserSession {
    private ws: WebSocket | null = null;
    private config: FileBrowserConfig | null = null;
    private callbacks: FileBrowserCallbacks;
    private fileItems: FileItem[] = [];
    private currentCommand: string = '';
    private commandBuffer: string = '';
    private isConnected: boolean = false;
    private isReady: boolean = false;

    constructor(callbacks: FileBrowserCallbacks = {}) {
        this.callbacks = callbacks;
    }

    connect(config: FileBrowserConfig) {
        if (this.ws?.readyState === WebSocket.OPEN) {
            this.ws.close();
        }

        this.config = config;
        const wsUrl = config.wsUrl || `${OMICS_FLOW_WS_URL}/script-runner`;
        this.ws = new WebSocket(wsUrl);

        this.ws.onopen = () => {
            this.isConnected = true;
            this.callbacks.onConnect?.();

            const message = {
                type: 'run',
                node_id: config.node_id,
                linux_user: config.linux_user,
                private_key: config.private_key,
                environment: config.environment || 'base',
                host: config.host || '127.0.0.1',
                port: config.port || 22,
            };

            this.ws!.send(JSON.stringify(message));
        };

        this.ws.onmessage = (event) => {
            const message = JSON.parse(event.data);

            switch (message.type) {
                case 'stdout':
                    this.handleStdout(message.data);
                    break;
                case 'stderr':
                    break;
                case 'status':
                    if (!this.isReady) {
                        this.isReady = true;
                        setTimeout(() => {
                            const initialPath = config.initialPath || `/nodes/${config.node_id}/`;
                            this.loadDirectory(initialPath, 0);
                        }, 300);
                    }
                    break;
                case 'error':
                    this.callbacks.onError?.(message.data);
                    break;
            }
        };

        this.ws.onerror = () => {
            this.callbacks.onError?.('WebSocket connection error');
        };

        this.ws.onclose = () => {
            this.isConnected = false;
            this.isReady = false;
            this.callbacks.onClose?.();
        };
    }

    private loadDirectory(path: string, depth: number) {
        if (!this.ws || !this.isReady) return;

        this.currentCommand = `ls_${path}_${depth}`;
        this.commandBuffer = '';

        const command = `ls -1p "${path}"\n`;
        this.ws.send(JSON.stringify({ type: 'input', data: command }));
    }

    private handleStdout(data: string) {
        this.commandBuffer += data;
        console.log("Received stdout:", data);

        if (data.includes('$') || data.includes('#') || data.includes('Singularity>')) {
            this.processLsOutput();
        }
    }

    private processLsOutput() {
        if (!this.currentCommand.startsWith('ls_')) return;

        const [, path, depthStr] = this.currentCommand.split('_');
        const depth = parseInt(depthStr, 10);


        const lines = this.commandBuffer
            .split(/\r?\n/)
            .map(line => line.replace(/\x1b\[[0-9;?]*[a-zA-Z]/g, '').trim())
            .filter(line => {
                if (!line) return false;
                if (line.startsWith('ls ')) return false;
                if (line.includes('$')) return false;
                if (line.includes('#')) return false;
                if (line.includes('Singularity>')) return false;
                if (line.match(/^[!@].*:/)) return false;
                return true;
            });

        console.log("filtered lines", lines);

        const newItems: FileItem[] = lines.map(line => {
            const isDirectory = line.endsWith('/');
            const name = isDirectory ? line.slice(0, -1) : line;
            const fullPath = path.endsWith('/') ? `${path}${name}` : `${path}/${name}`;

            return {
                name,
                fullPath,
                type: isDirectory ? 'directory' : 'file',
                isExpanded: false,
                depth
            };
        });

        if (depth === 0) {
            this.fileItems = newItems;
        } else {
            const parentPath = path.endsWith('/') ? path.slice(0, -1) : path;
            const parentIndex = this.fileItems.findIndex(item => item.fullPath === parentPath);
            if (parentIndex !== -1) {
                this.fileItems.splice(parentIndex + 1, 0, ...newItems);
            }
        }

        this.callbacks.onFileStructureChange?.([...this.fileItems]);
        this.currentCommand = '';
        this.commandBuffer = '';
    }

    toggleDirectory(item: FileItem) {
        const index = this.fileItems.findIndex(i => i.fullPath === item.fullPath);
        if (index === -1 || item.type !== 'directory') return;

        if (item.isExpanded) {
            this.fileItems[index].isExpanded = false;
            const childrenToRemove = this.fileItems.filter(
                i => i.depth > item.depth && i.fullPath.startsWith(item.fullPath + '/')
            );
            this.fileItems = this.fileItems.filter(i => !childrenToRemove.includes(i));
        } else {
            this.fileItems[index].isExpanded = true;
            const fullPath = item.fullPath.endsWith('/') ? item.fullPath : item.fullPath + '/';
            this.loadDirectory(fullPath, item.depth + 1);
        }

        this.callbacks.onFileStructureChange?.([...this.fileItems]);
    }

    getFileItems(): FileItem[] {
        return [...this.fileItems];
    }

    getIsConnected(): boolean {
        return this.isConnected;
    }

    close() {
        if (this.ws) {
            this.ws.close();
            this.ws = null;
        }
    }
}
