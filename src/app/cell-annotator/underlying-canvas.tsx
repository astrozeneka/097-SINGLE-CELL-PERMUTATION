"use client";

// ---------------------------------------------------------------------------
// ColorEncoder<T>
// Describes how to color each point entirely from the consumer side.
// The canvas auto-generates shader declarations and builds the interleaved
// buffer — the consumer only writes the domain-specific parts.
//
// colorGlsl has access to every declared attribute and uniform by name.
// The hue2rgb() helper is always available in colorGlsl.
//
// Example (color by cluster):
//   const enc: ColorEncoder<Cell> = {
//       attributes: [{ name: "a_cluster", size: 1, feed: d => d.cluster }],
//       uniforms:   [{ name: "u_n", type: "float", value: 14 }],
//       colorGlsl:  `return vec4(hue2rgb(a_cluster / u_n), 0.8);`,
//   };
// ---------------------------------------------------------------------------
export interface ColorEncoder<T> {
    // Extra per-point floats appended after x, y in the interleaved buffer.
    // Each attribute is auto-declared in the vertex shader and passed as a
    // varying to the fragment shader.
    attributes: {
        name: string;
        size: 1 | 2 | 3 | 4;
        feed: (d: T) => number | number[];
    }[];

    // Scalar/vector uniforms uploaded once before drawing.
    uniforms: {
        name: string;
        type: "float" | "vec2" | "vec3" | "vec4";
        value: number | number[];
    }[];

    // GLSL body of the color function — must contain a `return vec4(...)`.
    // May reference any declared attribute or uniform by name.
    // hue2rgb(h: float) -> vec3  is always available.
    colorGlsl: string;
}

// Zoom/pan state. Managed externally so an overlay canvas (e.g. drawing layer)
// can share the exact same transform without coupling to this component.
export interface Transform {
    x: number;
    y: number;
    scale: number;
}

interface UnderlyingCanvasParams<T> {
    data: T[];

    // Position accessors. The canvas auto-fits all points on init.
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;

    // Optional third dimension (e.g. t-SNE / UMAP z-axis).
    // Projected away for 2D display — reserved for future 3D support.
    zAccessor?: (d: T) => number;

    // Drives all coloring logic. Swap at runtime to recolor without
    // touching the canvas implementation (cluster → cell type → sample, etc.).
    colorEncoder: ColorEncoder<T>;

    // Controlled externally so sibling canvases stay in sync.
    transform: Transform;
}

export function UnderlyingCanvas<T>(params: UnderlyingCanvasParams<T>) {

}
