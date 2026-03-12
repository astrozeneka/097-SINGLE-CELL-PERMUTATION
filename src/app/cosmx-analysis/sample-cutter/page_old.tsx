"use client"
import { useEffect, useRef, useState } from "react"

// 12-color palette for clusters
const PALETTE: [number, number, number][] = [
    [0.94, 0.28, 0.28], [0.28, 0.61, 0.94], [0.28, 0.86, 0.56],
    [0.95, 0.66, 0.22], [0.71, 0.36, 0.90], [0.94, 0.43, 0.75],
    [0.40, 0.81, 0.81], [0.85, 0.85, 0.28], [0.56, 0.56, 0.56],
    [0.61, 0.36, 0.18], [0.94, 0.56, 0.28], [0.28, 0.94, 0.94],
]

const VERT = `
attribute vec2 a_pos;
attribute vec3 a_color;
uniform vec2 u_offset;
uniform float u_scale;
varying vec3 v_color;
void main() {
    vec2 p = (a_pos + u_offset) * u_scale * vec2(1.0, -1.0);
    gl_Position = vec4(p, 0.0, 1.0);
    gl_PointSize = 5.0;
    v_color = a_color;
}`

const FRAG = `
#extension GL_OES_standard_derivatives : enable
precision mediump float;
varying vec3 v_color;
void main() {
    float d = length(gl_PointCoord - vec2(0.5));
    float alpha = 1.0 - smoothstep(0.5 - fwidth(d), 0.5, d);
    gl_FragColor = vec4(v_color, alpha * 0.85);
}`

function makeShader(gl: WebGLRenderingContext, type: number, src: string) {
    const s = gl.createShader(type)!
    gl.shaderSource(s, src)
    gl.compileShader(s)
    return s
}

function parseCSV(csv: string): { x: number; y: number; cluster: string }[] {
    const lines = csv.trim().split("\n")
    const hdr = lines[0].split(",").map(h => h.trim())
    const xi = hdr.indexOf("x"), yi = hdr.indexOf("y"), ci = hdr.indexOf("cluster")
    return lines.slice(1).map(l => {
        const c = l.split(",")
        return { x: parseFloat(c[xi]), y: parseFloat(c[yi]), cluster: c[ci]?.trim() ?? "" }
    })
}

export default function SampleCutterPage() {
    const containerRef = useRef<HTMLDivElement>(null)
    const glCanvasRef = useRef<HTMLCanvasElement>(null)
    const overlayRef = useRef<HTMLCanvasElement>(null)
    const [geojson, setGeojson] = useState<object | null>(null)
    const [hasData, setHasData] = useState(false)

    // WebGL handles
    const glRef = useRef<WebGLRenderingContext | null>(null)
    const uOffRef = useRef<WebGLUniformLocation | null>(null)
    const uScaleRef = useRef<WebGLUniformLocation | null>(null)
    const cancelRafRef = useRef<() => void>(() => {})
    const dirtyRef = useRef(false)

    // View: clip_x = (data_x + ox) * s,  clip_y = -(data_y + oy) * s
    const view = useRef({ ox: 0, oy: 0, s: 1 })

    // Polygon vertices in data space
    const poly = useRef<{ x: number; y: number }[]>([])
    const closed = useRef(false)
    const drawOverlayRef = useRef<() => void>(() => {})

    useEffect(() => {
        const glCanvas = glCanvasRef.current!
        const overlay = overlayRef.current!
        const container = containerRef.current!

        const sizeCanvases = () => {
            glCanvas.width = container.clientWidth
            glCanvas.height = container.clientHeight
            overlay.width = glCanvas.width
            overlay.height = glCanvas.height
            dirtyRef.current = true
            drawOverlayRef.current()
        }
        sizeCanvases()
        const ro = new ResizeObserver(sizeCanvases)
        ro.observe(container)

        // Screen pixel → data coord
        const toData = (sx: number, sy: number) => {
            const { ox, oy, s } = view.current
            const cx = (sx / overlay.width) * 2 - 1
            const cy = 1 - (sy / overlay.height) * 2
            return { x: cx / s - ox, y: -cy / s - oy }
        }

        // Data coord → overlay screen pixel
        const toScreen = (dx: number, dy: number) => {
            const { ox, oy, s } = view.current
            const cx = (dx + ox) * s
            const cy = -(dy + oy) * s
            return { x: (cx + 1) / 2 * overlay.width, y: (-cy + 1) / 2 * overlay.height }
        }

        const drawOverlay = () => {
            const ctx = overlay.getContext("2d")!
            ctx.clearRect(0, 0, overlay.width, overlay.height)
            const pts = poly.current
            if (pts.length === 0) return
            ctx.strokeStyle = "rgba(255,255,255,0.9)"
            ctx.lineWidth = 2
            ctx.beginPath()
            const p0 = toScreen(pts[0].x, pts[0].y)
            ctx.moveTo(p0.x, p0.y)
            for (let i = 1; i < pts.length; i++) {
                const p = toScreen(pts[i].x, pts[i].y)
                ctx.lineTo(p.x, p.y)
            }
            if (closed.current) ctx.closePath()
            ctx.stroke()
            ctx.fillStyle = "white"
            for (const pt of pts) {
                const p = toScreen(pt.x, pt.y)
                ctx.beginPath()
                ctx.arc(p.x, p.y, 4, 0, Math.PI * 2)
                ctx.fill()
            }
        }
        drawOverlayRef.current = drawOverlay

        // Pan/zoom helpers
        const getPos = (e: MouseEvent) => {
            const rect = overlay.getBoundingClientRect()
            return {
                x: (e.clientX - rect.left) * (overlay.width / rect.width),
                y: (e.clientY - rect.top) * (overlay.height / rect.height),
            }
        }

        const updateGLView = () => {
            const gl = glRef.current
            if (!gl || !uOffRef.current || !uScaleRef.current) return
            gl.uniform2f(uOffRef.current, view.current.ox, view.current.oy)
            gl.uniform1f(uScaleRef.current, view.current.s)
            dirtyRef.current = true
        }

        const onWheel = (e: WheelEvent) => {
            e.preventDefault()
            const { x: sx, y: sy } = getPos(e)
            const f = e.deltaY < 0 ? 1.1 : 1 / 1.1
            const { ox, oy, s } = view.current
            const cx = (sx / overlay.width) * 2 - 1
            const cy = 1 - (sy / overlay.height) * 2
            // Keep the data point under cursor fixed after zoom
            const dataX = cx / s - ox, dataY = -cy / s - oy
            const ns = s * f
            view.current = { ox: cx / ns - dataX, oy: -cy / ns - dataY, s: ns }
            updateGLView()
            drawOverlay()
        }

        let dragging = false, dragPrev = { x: 0, y: 0 }, dragMoved = false

        const onMouseDown = (e: MouseEvent) => {
            dragging = true
            dragMoved = false
            dragPrev = getPos(e)
        }

        const onMouseMove = (e: MouseEvent) => {
            if (!dragging) return
            const pos = getPos(e)
            const dx = pos.x - dragPrev.x, dy = pos.y - dragPrev.y
            if (Math.abs(dx) > 3 || Math.abs(dy) > 3) dragMoved = true
            if (!dragMoved) return
            // Pan: offset shifts so the point under cursor follows the mouse
            const { s } = view.current
            view.current.ox += dx / overlay.width * 2 / s
            view.current.oy += dy / overlay.height * 2 / s
            dragPrev = pos
            updateGLView()
            drawOverlay()
        }

        const onMouseUp = (e: MouseEvent) => {
            if (dragging && !dragMoved && !closed.current) {
                const pos = getPos(e)
                poly.current.push(toData(pos.x, pos.y))
                drawOverlay()
            }
            dragging = false
        }

        const onKeyDown = (e: KeyboardEvent) => {
            if (e.key.toLowerCase() !== "c" || closed.current || poly.current.length < 3) return
            closed.current = true
            const ring = [...poly.current, poly.current[0]].map(p => [p.x, p.y])
            setGeojson({ type: "Feature", geometry: { type: "Polygon", coordinates: [ring] }, properties: {} })
            drawOverlay()
        }

        overlay.addEventListener("wheel", onWheel, { passive: false })
        overlay.addEventListener("mousedown", onMouseDown)
        window.addEventListener("mousemove", onMouseMove)
        window.addEventListener("mouseup", onMouseUp)
        window.addEventListener("keydown", onKeyDown)

        return () => {
            cancelRafRef.current()
            ro.disconnect()
            overlay.removeEventListener("wheel", onWheel)
            overlay.removeEventListener("mousedown", onMouseDown)
            window.removeEventListener("mousemove", onMouseMove)
            window.removeEventListener("mouseup", onMouseUp)
            window.removeEventListener("keydown", onKeyDown)
        }
    }, [])

    const initWebGL = (data: { x: number; y: number; cluster: string }[]) => {
        const glCanvas = glCanvasRef.current!
        const gl = glCanvas.getContext("webgl", { antialias: true, alpha: false })!
        gl.getExtension("OES_standard_derivatives")
        glRef.current = gl

        const prog = gl.createProgram()!
        gl.attachShader(prog, makeShader(gl, gl.VERTEX_SHADER, VERT))
        gl.attachShader(prog, makeShader(gl, gl.FRAGMENT_SHADER, FRAG))
        gl.linkProgram(prog)
        gl.useProgram(prog)

        const n = data.length
        const pos = new Float32Array(n * 2)
        const colors = new Float32Array(n * 3)
        const clusterIdx = new Map<string, number>()
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity

        for (let i = 0; i < n; i++) {
            const { x, y, cluster } = data[i]
            pos[i * 2] = x; pos[i * 2 + 1] = y
            if (!clusterIdx.has(cluster)) clusterIdx.set(cluster, clusterIdx.size)
            const [r, g, b] = PALETTE[clusterIdx.get(cluster)! % PALETTE.length]
            colors[i * 3] = r; colors[i * 3 + 1] = g; colors[i * 3 + 2] = b
            if (x < minX) minX = x; if (x > maxX) maxX = x
            if (y < minY) minY = y; if (y > maxY) maxY = y
        }

        const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2
        const initScale = 1.8 / Math.max(maxX - minX, maxY - minY)
        view.current = { ox: -cx, oy: -cy, s: initScale }

        const mkBuf = (arr: Float32Array, attr: string, size: number) => {
            const buf = gl.createBuffer()!
            gl.bindBuffer(gl.ARRAY_BUFFER, buf)
            gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW)
            const loc = gl.getAttribLocation(prog, attr)
            gl.enableVertexAttribArray(loc)
            gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0)
        }
        mkBuf(pos, "a_pos", 2)
        mkBuf(colors, "a_color", 3)

        uOffRef.current = gl.getUniformLocation(prog, "u_offset")
        uScaleRef.current = gl.getUniformLocation(prog, "u_scale")
        gl.uniform2f(uOffRef.current, -cx, -cy)
        gl.uniform1f(uScaleRef.current, initScale)
        gl.enable(gl.BLEND)
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
        gl.clearColor(0.05, 0.05, 0.08, 1)

        dirtyRef.current = true
        let rafId: number
        const render = () => {
            if (dirtyRef.current) {
                dirtyRef.current = false
                gl.viewport(0, 0, glCanvas.width, glCanvas.height)
                gl.clear(gl.COLOR_BUFFER_BIT)
                gl.drawArrays(gl.POINTS, 0, n)
            }
            rafId = requestAnimationFrame(render)
        }
        rafId = requestAnimationFrame(render)
        cancelRafRef.current = () => cancelAnimationFrame(rafId)
        setHasData(true)
    }

    const handleFile = (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0]
        if (!file) return
        const reader = new FileReader()
        reader.onload = ev => initWebGL(parseCSV(ev.target!.result as string))
        reader.readAsText(file)
    }

    return (
        <div ref={containerRef} style={{ width: "100vw", height: "100vh", position: "relative", background: "#0d0d12" }}>
            <canvas ref={glCanvasRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%" }} />
            <canvas ref={overlayRef} style={{ position: "absolute", inset: 0, width: "100%", height: "100%", cursor: "crosshair" }} />

            {!hasData && (
                <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center" }}>
                    <label style={{ padding: "12px 24px", background: "#1e1e2e", color: "white", borderRadius: 6, cursor: "pointer", border: "1px solid #444" }}>
                        Load CSV (x, y, cluster)
                        <input type="file" accept=".csv" onChange={handleFile} style={{ display: "none" }} />
                    </label>
                </div>
            )}

            {hasData && !geojson && (
                <div style={{ position: "absolute", top: 12, left: 12, color: "rgba(255,255,255,0.45)", fontSize: 12, pointerEvents: "none" }}>
                    Click to place points · Press C to close polygon (min 3 points)
                </div>
            )}

            {geojson && (
                <div style={{ position: "absolute", bottom: 16, left: 16, right: 16, maxHeight: "35vh", overflow: "auto", background: "rgba(0,0,0,0.85)", borderRadius: 6, padding: 12, border: "1px solid #444" }}>
                    <pre style={{ color: "#7cf", margin: 0, fontSize: 11, userSelect: "all" }}>
                        {JSON.stringify(geojson, null, 2)}
                    </pre>
                </div>
            )}
        </div>
    )
}
