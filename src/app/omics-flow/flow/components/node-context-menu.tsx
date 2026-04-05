

export function NodeContextMenu({ contextMenu, onAction }: { contextMenu: { x: number; y: number } | null; onAction: (action: 'copy' | 'cut' | 'delete') => void }) {
    if (!contextMenu) return null;

    return (
        <div
            style={{
                position: "fixed",
                left: `${contextMenu.x}px`,
                top: `${contextMenu.y}px`,
                background: "white",
                border: "1px solid #ccc",
                borderRadius: "4px",
                boxShadow: "0 2px 8px rgba(0,0,0,0.15)",
                zIndex: 10000,
                minWidth: "150px"
            }}
        >
            <div
                onClick={() => onAction('copy')}
                style={{
                    padding: "8px 16px",
                    cursor: "pointer",
                    borderBottom: "1px solid #eee"
                }}
                onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                onMouseLeave={(e) => e.currentTarget.style.background = "white"}
            >
                Copy clone
            </div>
            <div
                onClick={() => onAction('cut')}
                style={{
                    padding: "8px 16px",
                    cursor: "pointer",
                    borderBottom: "1px solid #eee"
                }}
                onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                onMouseLeave={(e) => e.currentTarget.style.background = "white"}
            >
                Cut
            </div>
            <div
                onClick={() => onAction('delete')}
                style={{
                    padding: "8px 16px",
                    cursor: "pointer",
                    color: "#d32f2f"
                }}
                onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                onMouseLeave={(e) => e.currentTarget.style.background = "white"}
            >
                Delete
            </div>
        </div>
    );
}