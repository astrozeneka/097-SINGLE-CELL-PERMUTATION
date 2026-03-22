import fs from "fs";
import path from "path";
import { SampleCutterClient } from "./sample-cutter-client";

function loadData(file: string): { data: number[]; numClusters: number; dataBounds: { minX: number; maxX: number; minY: number; maxY: number } } {
    const csv = fs.readFileSync(path.join(process.cwd(), "public/sample-datas/spatial-cosmx", file), "utf8");
    const data: number[] = [];
    for (const line of csv.split("\n").slice(1)) {
        if (!line) continue;
        const parts = line.split(",");
        data.push(parseFloat(parts[1]), parseFloat(parts[2]), parseFloat(parts[3]));
    }

    // Invert x coordinates
    let maxX = -Infinity;
    for (let i = 0; i < data.length; i += 3) maxX = Math.max(maxX, data[i]);
    for (let i = 0; i < data.length; i += 3) data[i] = maxX - data[i];

    // Count distinct cluster values
    const clusters = new Set<number>();
    for (let i = 2; i < data.length; i += 3) clusters.add(data[i]);

    let minX = Infinity, minY = Infinity, maxY = -Infinity;
    maxX = -Infinity;
    for (let i = 0; i < data.length; i += 3) {
        if (data[i] < minX) minX = data[i]; if (data[i] > maxX) maxX = data[i];
        if (data[i + 1] < minY) minY = data[i + 1]; if (data[i + 1] > maxY) maxY = data[i + 1];
    }

    return { data, numClusters: clusters.size, dataBounds: { minX, maxX, minY, maxY } };
}

function loadPolygons(sampleName: string, maxX: number) {
    const dir = path.join(process.cwd(), "public/sample-datas/spatial-cosmx/geojson", sampleName);
    if (!fs.existsSync(dir)) return [];
    return fs.readdirSync(dir)
        .filter(f => f.endsWith(".geojson"))
        .map(f => {
            const { geometry } = JSON.parse(fs.readFileSync(path.join(dir, f), "utf8"));
            const ring: number[][] = geometry.coordinates[0];
            return { name: f.slice(0, -8), points: ring.slice(0, -1).map(([x, y]) => ({ x: maxX - x, y })) };
        });
}

export default async function SampleCutterPage({ searchParams }: { searchParams: Promise<{ sample?: string }> }) {
    const { sample } = await searchParams;
    const sampleName = sample === "LARC_A" ? "LARC_A" : "LARC_B";
    const { data, numClusters, dataBounds } = loadData(`${sampleName}.csv`);
    const savedPolygons = loadPolygons(sampleName, dataBounds.maxX);
    return <SampleCutterClient data={data} numClusters={numClusters} savedPolygons={savedPolygons} dataBounds={dataBounds} />;
}
