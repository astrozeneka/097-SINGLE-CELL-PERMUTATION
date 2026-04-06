import { useState, useEffect } from "react";

interface CreateNodeDialogProps {
    isOpen: boolean;
    onClose: () => void;
    onCreate: (name: string, description: string) => void;
}

export function CreateNodeDialog({ isOpen, onClose, onCreate }: CreateNodeDialogProps) {
    const [name, setName] = useState("");
    const [description, setDescription] = useState("");
    const [dirname, setDirname] = useState("");

    useEffect(() => {
        if (name) {
            const sanitized = name.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
            const randomSuffix = Math.floor(Math.random() * 100000).toString().padStart(5, '0');
            setDirname(`${sanitized}-${randomSuffix}`);
        } else {
            setDirname("");
        }
    }, [name]);

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (name.trim()) {
            onCreate(name.trim(), description.trim());
            setName("");
            setDescription("");
            setDirname("");
            onClose();
        }
    };

    const handleCancel = () => {
        setName("");
        setDescription("");
        setDirname("");
        onClose();
    };

    if (!isOpen) return null;

    return (
        <div
            style={{
                position: "fixed",
                top: 0,
                left: 0,
                right: 0,
                bottom: 0,
                background: "rgba(0, 0, 0, 0.5)",
                display: "flex",
                alignItems: "center",
                justifyContent: "center",
                zIndex: 10001
            }}
            onClick={handleCancel}
        >
            <div
                style={{
                    background: "white",
                    borderRadius: "8px",
                    padding: "24px",
                    minWidth: "400px",
                    maxWidth: "500px",
                    boxShadow: "0 4px 16px rgba(0, 0, 0, 0.2)"
                }}
                onClick={(e) => e.stopPropagation()}
            >
                <h2 style={{ margin: "0 0 20px 0", fontSize: "20px" }}>Create New Node</h2>
                <form onSubmit={handleSubmit}>
                    <div style={{ marginBottom: "16px" }}>
                        <label
                            htmlFor="node-name"
                            style={{
                                display: "block",
                                marginBottom: "6px",
                                fontSize: "14px",
                                fontWeight: "500"
                            }}
                        >
                            Name <span style={{ color: "red" }}>*</span>
                        </label>
                        <input
                            id="node-name"
                            type="text"
                            value={name}
                            onChange={(e) => setName(e.target.value)}
                            placeholder="Enter node name"
                            style={{
                                width: "100%",
                                padding: "8px 12px",
                                fontSize: "14px",
                                border: "1px solid #ccc",
                                borderRadius: "4px",
                                boxSizing: "border-box"
                            }}
                            required
                            autoFocus
                        />
                    </div>

                    <div style={{ marginBottom: "16px" }}>
                        <label
                            htmlFor="node-dirname"
                            style={{
                                display: "block",
                                marginBottom: "6px",
                                fontSize: "14px",
                                fontWeight: "500"
                            }}
                        >
                            Directory Name (auto-generated)
                        </label>
                        <input
                            id="node-dirname"
                            type="text"
                            value={dirname}
                            readOnly
                            style={{
                                width: "100%",
                                padding: "8px 12px",
                                fontSize: "14px",
                                border: "1px solid #e0e0e0",
                                borderRadius: "4px",
                                boxSizing: "border-box",
                                background: "#f5f5f5",
                                color: "#666",
                                cursor: "not-allowed"
                            }}
                        />
                    </div>

                    <div style={{ marginBottom: "20px" }}>
                        <label
                            htmlFor="node-description"
                            style={{
                                display: "block",
                                marginBottom: "6px",
                                fontSize: "14px",
                                fontWeight: "500"
                            }}
                        >
                            Description
                        </label>
                        <textarea
                            id="node-description"
                            value={description}
                            onChange={(e) => setDescription(e.target.value)}
                            placeholder="Enter node description (optional)"
                            rows={3}
                            style={{
                                width: "100%",
                                padding: "8px 12px",
                                fontSize: "14px",
                                border: "1px solid #ccc",
                                borderRadius: "4px",
                                boxSizing: "border-box",
                                resize: "vertical",
                                fontFamily: "inherit"
                            }}
                        />
                    </div>

                    <div style={{ display: "flex", justifyContent: "flex-end", gap: "12px" }}>
                        <button
                            type="button"
                            onClick={handleCancel}
                            style={{
                                padding: "8px 16px",
                                fontSize: "14px",
                                border: "1px solid #ccc",
                                borderRadius: "4px",
                                background: "white",
                                cursor: "pointer"
                            }}
                        >
                            Cancel
                        </button>
                        <button
                            type="submit"
                            disabled={!name.trim()}
                            style={{
                                padding: "8px 16px",
                                fontSize: "14px",
                                border: "none",
                                borderRadius: "4px",
                                background: name.trim() ? "#4A90E2" : "#ccc",
                                color: "white",
                                cursor: name.trim() ? "pointer" : "not-allowed"
                            }}
                        >
                            Create Node
                        </button>
                    </div>
                </form>
            </div>
        </div>
    );
}
