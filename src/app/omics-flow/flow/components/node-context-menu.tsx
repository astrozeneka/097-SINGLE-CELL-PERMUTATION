

export function NodeContextMenu({
    contextMenu,
    onAction
}: {
    contextMenu: { x: number; y: number; type: 'node' | 'canvas'; nodes: any[] } | null;
    onAction: (action: 'copy' | 'cut' | 'delete' | 'paste', cursorPosition?: { x: number; y: number }) => void;
}) {
    if (!contextMenu) return null;

    const isCanvasMenu = contextMenu.type === 'canvas';

    const handleAction = (action: 'copy' | 'cut' | 'delete' | 'paste', cursorPosition?: { x: number; y: number }) => {
        return (e: React.MouseEvent) => {
            e.stopPropagation();
            onAction(action, cursorPosition);
        };
    };

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
            onClick={(e) => e.stopPropagation()}
            onMouseDown={(e) => e.stopPropagation()}
            onMouseUp={(e) => e.stopPropagation()}
        >
            {isCanvasMenu ? (
                // Canvas menu: only show Paste
                <div
                    onClick={handleAction('paste', { x: contextMenu.x, y: contextMenu.y })}
                    style={{
                        padding: "8px 16px",
                        cursor: "pointer"
                    }}
                    onMouseEnter={(e) => e.currentTarget.style.background = "#f5f5f5"}
                    onMouseLeave={(e) => e.currentTarget.style.background = "white"}
                >
                    Paste
                </div>
            ) : (
                // Node menu: show Copy, Cut, Delete
                <>
                    <div
                        onClick={handleAction('copy')}
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
                        onClick={handleAction('cut')}
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
                        onClick={handleAction('delete')}
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
                </>
            )}
        </div>
    );
}