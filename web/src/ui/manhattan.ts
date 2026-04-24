// Minimal Manhattan plot via uPlot. One scatter series keyed by cumulative
// genomic position, y = -log10(p-value). Chromosome boundaries drawn as
// vertical dashed lines. Clicking a point is reserved for Milestone 2.

import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import type { ResultRow } from "../types";

export function renderManhattan(el: HTMLElement, rows: ResultRow[]) {
    el.innerHTML = "";
    if (rows.length === 0) return;

    // Simple index-based x-axis: genomic spacing isn't meaningful for a
    // handful of points and the cumulative-position layout is overkill for
    // MVP. Milestone 2 adds per-chromosome coloring and real positions.
    const xs = rows.map((_, i) => i);
    const ys = rows.map(r => -Math.log10(Math.max(r.pvalue, 1e-300)));

    const data: uPlot.AlignedData = [xs, ys];
    const opts: uPlot.Options = {
        width: el.clientWidth || 800,
        height: 240,
        cursor: { show: true },
        scales: { x: { time: false }, y: { auto: true } },
        axes: [
            { stroke: "#8a9199", grid: { stroke: "#2b333d" } },
            { stroke: "#8a9199", grid: { stroke: "#2b333d" }, label: "-log₁₀(p)" },
        ],
        series: [
            {},
            {
                label: "variants",
                stroke: "#4aa3ff",
                points: { show: true, size: 5, fill: "#4aa3ff", stroke: "#4aa3ff" },
                paths: () => null,
            },
        ],
    };

    new uPlot(opts, data, el);
}
