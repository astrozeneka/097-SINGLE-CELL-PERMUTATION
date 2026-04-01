import { ColorEncoder, Transform } from "../underlying-canvas";

export function Scatter3DV2<T>({ data, xAccessor, yAccessor, zAccessor, colorEncoder, size, transform, setTransform, onReady }: {
    data: T[];
    xAccessor: (d: T) => number;
    yAccessor: (d: T) => number;
    zAccessor: (d: T) => number;
    colorEncoder: ColorEncoder<T>;
    size: { w: number; h: number };
    transform: Transform;
    setTransform: (transform: Transform) => void;
    onReady: () => void;
}) {
    return <div>hello 3d</div>;
}