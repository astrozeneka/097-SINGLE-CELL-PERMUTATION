

export default function SubsetSelectorV2_1({ subsets, subset, onSubsetChange }: { subsets: string[]; subset: string | null; onSubsetChange: (s: string) => void }) {
    return (
        <div className="thin-scrollbar" style={{
            position: "absolute", bottom: 0, left: "50%", transform: "translateX(-50%)",
            display: "flex", gap: 6, padding: "6px 10px",
            background: "rgba(0,0,0,0.55)", backdropFilter: "blur(6px)",
            borderRadius: 10, overflowX: "auto", maxWidth: "100%",
            boxShadow: "0 2px 12px rgba(0,0,0,0.4)",
        }}>
            {subsets.map(p => (
                <button key={p} onClick={() => onSubsetChange(p)} style={{
                    flexShrink: 0, padding: "4px 12px", borderRadius: 6, border: "none",
                    cursor: "pointer", fontSize: 12, fontFamily: "monospace",
                    background: p === subset ? "rgba(255,255,255,0.9)" : "rgba(255,255,255,0.12)",
                    color: p === subset ? "#111" : "#ddd",
                    fontWeight: p === subset ? 700 : 400,
                }}>
                    {p}
                </button>
            ))}
        </div>
    );
}