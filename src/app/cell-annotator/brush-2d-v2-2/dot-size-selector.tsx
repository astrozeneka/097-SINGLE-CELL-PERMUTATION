export default function DotSizeSelector({ dotSize, onDotSizeChange }: { dotSize: number; onDotSizeChange: (s: number) => void }) {
    return (
        <div style={{
            position: "absolute", top: 12, left: 12,
            display: "flex", alignItems: "center", gap: 8, padding: "6px 10px",
            background: "rgba(0,0,0,0.55)", backdropFilter: "blur(6px)",
            borderRadius: 10, boxShadow: "0 2px 12px rgba(0,0,0,0.4)",
        }}>
            <span style={{ fontSize: 11, fontFamily: "monospace", color: "#aaa" }}>dot</span>
            <button onClick={() => onDotSizeChange(Math.max(1, dotSize - 1))} style={btnStyle}>−</button>
            <span style={{ fontSize: 12, fontFamily: "monospace", color: "#eee", minWidth: 20, textAlign: "center" }}>{dotSize}</span>
            <button onClick={() => onDotSizeChange(dotSize + 1)} style={btnStyle}>+</button>
        </div>
    );
}

const btnStyle: React.CSSProperties = {
    width: 24, height: 24, borderRadius: 6, border: "none",
    background: "rgba(255,255,255,0.12)", color: "#ddd",
    cursor: "pointer", fontSize: 16, lineHeight: 1,
    display: "flex", alignItems: "center", justifyContent: "center",
};
