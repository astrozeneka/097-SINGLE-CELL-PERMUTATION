

// Zoom/pan state. Managed externally so an overlay canvas (e.g. drawing layer)
// can share the exact same transform without coupling to this component.
export interface Transform {
    x: number;
    y: number;
    scale: number;
}