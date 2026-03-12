import fs from "fs";
import path from "path";
import { SampleCutterClient } from "./sample-cutter-client";

function loadData(file: string): number[] {
    const csv = fs.readFileSync(path.join(process.cwd(), "public/sample-datas/spatial-cosmx", file), "utf8");
    const data: number[] = [];
    for (const line of csv.split("\n").slice(1)) {
        if (!line) continue;
        const parts = line.split(",");
        data.push(parseFloat(parts[1]), parseFloat(parts[2]), parseFloat(parts[3]));
    }
    return data;
}

export default function SampleCutterPage() {
    const data = loadData("LARC_B.csv");
    return <SampleCutterClient data={data} />;
}
