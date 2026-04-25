// Manhattan plot via uPlot. Real genomic positions on x (cumulative
// across chromosomes), -log10(p) on y, alternating colors per chrom,
// hover tooltip, drag-to-zoom (uPlot built-in), and click-to-drill-down.

import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import type { ResultRow } from "../types";

const CHROM_COLORS = ["#3c7ccc", "#a62a2a"];   // alternates across chroms
const SUGG_LINE   = "#3f7a4a";                  // suggestive threshold
const SIG_LINE    = "#a62a2a";                  // genome-wide threshold

// Sort key so 1..22, X, Y, MT come out in order.
function chromOrder(c: string): number {
    const stripped = c.replace(/^chr/i, "");
    if (/^\d+$/.test(stripped)) return parseInt(stripped, 10);
    if (stripped.toUpperCase() === "X")  return 23;
    if (stripped.toUpperCase() === "Y")  return 24;
    if (stripped.toUpperCase() === "MT" || stripped.toUpperCase() === "M") return 25;
    return 100 + stripped.charCodeAt(0);
}

interface ChromBand {
    chrom: string;
    start: number;        // cumulative x at start of chromosome
    end: number;          // cumulative x at end
    midpoint: number;     // for axis label
}

export function renderManhattan(
    el: HTMLElement,
    rows: ResultRow[],
    onPointClick?: (rowIdx: number, row: ResultRow) => void,
) {
    el.innerHTML = "";
    if (rows.length === 0) return;

    // Group rows by chromosome (preserve original row index for callbacks).
    const byChrom = new Map<string, { row: ResultRow; idx: number }[]>();
    rows.forEach((r, i) => {
        let arr = byChrom.get(r.chrom);
        if (!arr) { arr = []; byChrom.set(r.chrom, arr); }
        arr.push({ row: r, idx: i });
    });
    const chroms = Array.from(byChrom.keys()).sort((a, b) => chromOrder(a) - chromOrder(b));

    // Compute cumulative offsets so points across chromosomes don't overlap.
    // Use observed pos range per chromosome plus a small gap.
    const GAP = 1e6;          // 1Mb visual gap between chromosomes
    const bands: ChromBand[] = [];
    let cursor = 0;
    for (const c of chroms) {
        const arr = byChrom.get(c)!;
        let lo = Infinity, hi = -Infinity;
        for (const x of arr) { lo = Math.min(lo, x.row.pos); hi = Math.max(hi, x.row.pos); }
        if (lo === Infinity) { lo = 0; hi = 0; }
        const span = Math.max(1, hi - lo);
        const start = cursor;
        const end = cursor + span;
        bands.push({ chrom: c, start, end, midpoint: (start + end) / 2 });
        cursor = end + GAP;
    }
    const totalSpan = cursor - GAP;

    // One series per chromosome so each gets its own color. Each series's
    // y array is full-length with NaN where the row isn't in that chrom.
    // For small/medium VCFs (<<100k rows) this is plenty fast.
    const xs: number[] = new Array(rows.length);
    const ysByChrom: number[][] = chroms.map(() => new Array(rows.length).fill(NaN));
    chroms.forEach((c, ci) => {
        const arr = byChrom.get(c)!;
        const lo = Math.min(...arr.map(a => a.row.pos));
        for (const { row, idx } of arr) {
            xs[idx] = bands[ci].start + (row.pos - lo);
            ysByChrom[ci][idx] = -Math.log10(Math.max(row.pvalue, 1e-300));
        }
    });

    const data: uPlot.AlignedData = [xs, ...ysByChrom];

    // Tooltip element overlaid on the plot.
    const tooltip = document.createElement("div");
    tooltip.style.cssText = `
        position:absolute; pointer-events:none; z-index:5;
        background:#fdfbf5; border:1px solid #b9ab8d; border-radius:4px;
        padding:4px 8px; font:11px/1.3 ui-monospace, monospace;
        color:#2d0600; box-shadow:0 2px 6px rgba(0,0,0,0.08);
        display:none; white-space:nowrap;`;
    el.style.position = "relative";
    el.appendChild(tooltip);

    const opts: uPlot.Options = {
        width: el.clientWidth || 800,
        height: 240,
        cursor: {
            show: true,
            // Drag-to-zoom is uPlot default; just ensure x-axis drag is on.
            drag: { x: true, y: false, setScale: true },
        },
        scales: {
            x: { time: false, range: () => [0, totalSpan] },
            y: { auto: true },
        },
        axes: [
            {
                stroke: "#6f6354",
                // Custom tick values: use band midpoints, label with chrom name.
                values: (_self, _splits) => bands.map(b => b.chrom.replace(/^chr/i, "")),
                splits: () => bands.map(b => b.midpoint),
                grid: { show: false },
            },
            { stroke: "#6f6354", label: "-log₁₀(p)" },
        ],
        series: [
            {},
            ...chroms.map((c, ci) => ({
                label: c,
                stroke: CHROM_COLORS[ci % CHROM_COLORS.length],
                points: {
                    show: true,
                    size: 5,
                    fill: CHROM_COLORS[ci % CHROM_COLORS.length],
                    stroke: CHROM_COLORS[ci % CHROM_COLORS.length],
                },
                paths: () => null,
            })),
        ],
        hooks: {
            setCursor: [(u) => {
                const idx = u.cursor.idx;
                if (idx == null || idx < 0 || idx >= rows.length) {
                    tooltip.style.display = "none";
                    return;
                }
                const r = rows[idx];
                tooltip.innerHTML =
                    `<strong>${r.chrom}:${r.pos}</strong> ${r.ref}&gt;${r.alt}` +
                    `<br>p = ${r.pvalue.toExponential(3)}` +
                    `<br><span style="color:#6f6354">click for sample detail</span>`;
                const left = u.cursor.left ?? 0;
                const top  = u.cursor.top  ?? 0;
                tooltip.style.display = "block";
                tooltip.style.left = `${left + 12}px`;
                tooltip.style.top  = `${top + 12}px`;
            }],
        },
        plugins: [
            // Suggestive (1e-5) and genome-wide (5e-8) threshold lines.
            {
                hooks: {
                    draw: [(u) => {
                        const ctx = u.ctx;
                        const drawLine = (yVal: number, color: string, dash: number[]) => {
                            const yPx = u.valToPos(yVal, "y", true);
                            const x0 = u.bbox.left;
                            const x1 = u.bbox.left + u.bbox.width;
                            ctx.save();
                            ctx.strokeStyle = color;
                            ctx.setLineDash(dash);
                            ctx.beginPath();
                            ctx.moveTo(x0, yPx); ctx.lineTo(x1, yPx);
                            ctx.stroke();
                            ctx.restore();
                        };
                        const yMax = u.scales.y.max;
                        if (yMax != null && yMax >= 5)  drawLine(5,         SUGG_LINE, [4, 4]);
                        if (yMax != null && yMax >= 7.3) drawLine(7.30103,  SIG_LINE,  [4, 4]);
                    }],
                },
            },
        ],
    };

    const u = new uPlot(opts, data, el);

    // Resize on container changes (e.g. window resize, tab switch).
    const ro = new ResizeObserver(() => {
        u.setSize({ width: el.clientWidth, height: el.clientHeight || 240 });
    });
    ro.observe(el);

    if (onPointClick) {
        const over = el.querySelector<HTMLElement>(".u-over");
        if (over) {
            over.addEventListener("click", () => {
                const idx = u.cursor.idx;
                if (idx == null || idx < 0 || idx >= rows.length) return;
                onPointClick(idx, rows[idx]);
            });
        }
    }
}
