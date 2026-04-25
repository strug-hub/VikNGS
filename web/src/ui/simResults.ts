// Simulation result rendering: power curve, QQ plot, histogram.
// Mirrors SimPlotWindowPlotter.cpp (Qt). Plots fill available vertical
// space and re-fit on container resize.

import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import type { SimGroup, SimResultRow, SimRunRequest } from "../types";

const COLORS = ["#3c7ccc", "#a62a2a", "#3f7a4a", "#b68a00", "#8c46c6", "#2d8ba3"];

export interface SimRenderInput {
    rows: SimResultRow[];
    steps: number;
    family: "binomial" | "normal";
    /** Original group definitions, used by the sample-info side table. */
    groups: SimGroup[];
    summary: { variants: number; processingTime: number; evaluationTime: number };
}

function cohortLabel(g: SimGroup, family: "binomial" | "normal"): string {
    if (family === "normal") return "normal";
    return g.isCase ? "case" : "control";
}

function fmtDepth(g: SimGroup): string {
    return `${g.meanDepth.toFixed(1)} ± ${g.sdDepth.toFixed(1)} (${g.readDepth})`;
}

function sampleSizeAtStep(g: SimGroup, step: number): number {
    return g.n + g.nIncrement * step;
}

function rangeText(g: SimGroup, steps: number): string {
    if (g.nIncrement === 0 || steps <= 1) return String(g.n);
    return `${g.n} – ${sampleSizeAtStep(g, steps - 1)}`;
}

export function renderSampleTable(
    host: HTMLElement,
    groups: SimGroup[],
    family: "binomial" | "normal",
    step: number | null,    // null = show range across all steps
    steps: number,
) {
    const stepLabel = step === null
        ? (steps > 1 ? "min – max" : "")
        : `step ${step + 1}`;
    let html = `<h3>Sample info${stepLabel ? ` (${stepLabel})` : ""}</h3><table><thead><tr>
        <th>n</th><th>cohort</th><th>depth</th><th>err</th>
    </tr></thead><tbody>`;
    for (const g of groups) {
        const n = step === null ? rangeText(g, steps) : sampleSizeAtStep(g, step);
        html += `<tr>
            <td>${n}</td>
            <td>${cohortLabel(g, family)}</td>
            <td>${fmtDepth(g)}</td>
            <td>${g.errorRate.toFixed(3)}</td>
        </tr>`;
    }
    html += "</tbody></table>";
    host.innerHTML = html;
}

// Active uPlot instances and their ResizeObservers. Tracked so we can
// destroy them cleanly on re-render / tab switch.
const activePlots: { host: HTMLElement; plot: uPlot; ro: ResizeObserver }[] = [];

function destroyAllPlots() {
    for (const { plot, ro } of activePlots) { plot.destroy(); ro.disconnect(); }
    activePlots.length = 0;
}

function mountPlot(host: HTMLElement, opts: uPlot.Options, data: uPlot.AlignedData) {
    // Drop any plot previously attached to this exact host.
    for (let i = activePlots.length - 1; i >= 0; i--) {
        if (activePlots[i].host === host) {
            activePlots[i].plot.destroy();
            activePlots[i].ro.disconnect();
            activePlots.splice(i, 1);
        }
    }
    host.innerHTML = "";

    const measure = () => ({
        width: Math.max(50, host.clientWidth),
        height: Math.max(50, host.clientHeight),
    });
    const { width, height } = measure();
    const plot = new uPlot({ ...opts, width, height }, data, host);
    const ro = new ResizeObserver(() => plot.setSize(measure()));
    ro.observe(host);
    activePlots.push({ host, plot, ro });
}

export function renderSimResults(
    summaryEl: HTMLElement,
    powerEl: HTMLElement,
    qqEl: HTMLElement,
    histEl: HTMLElement,
    sampleEl: HTMLElement,
    input: SimRenderInput,
) {
    const { rows, steps, summary, groups, family } = input;
    summaryEl.innerHTML = `
      <strong>${summary.variants}</strong> variants &middot;
      processing <strong>${summary.processingTime.toFixed(2)}s</strong> &middot;
      evaluation <strong>${summary.evaluationTime.toFixed(2)}s</strong> &middot;
      <strong>${rows.length}</strong> p-values
    `;
    clearSimResults(summaryEl, powerEl, qqEl, histEl, sampleEl, /*keepSummary=*/true);

    // Sample-info side panel: defaults to the range across all steps;
    // hover updates to the specific step under the cursor.
    renderSampleTable(sampleEl, groups, family, steps > 1 ? null : 0, steps);

    // Group rows by (statName, genotypeSource) → list of { sampleSize, p }.
    const byKey = new Map<string, { statName: string; genotypeSource: string; pointsByStep: Map<number, number[]>; sampleSizeByStep: Map<number, number> }>();
    for (const r of rows) {
        const key = `${r.statName}/${r.genotypeSource}`;
        let entry = byKey.get(key);
        if (!entry) {
            entry = { statName: r.statName, genotypeSource: r.genotypeSource,
                      pointsByStep: new Map(), sampleSizeByStep: new Map() };
            byKey.set(key, entry);
        }
        let arr = entry.pointsByStep.get(r.stepIdx);
        if (!arr) { arr = []; entry.pointsByStep.set(r.stepIdx, arr); }
        arr.push(r.pvalue);
        entry.sampleSizeByStep.set(r.stepIdx, r.sampleSize);
    }

    const keys = Array.from(byKey.keys()).sort();
    const ALPHA = 0.05;

    function computePower(ps: number[], alpha: number) {
        if (ps.length === 0) return 0;
        let k = 0;
        for (const p of ps) if (p < alpha) k++;
        return k / ps.length;
    }

    // Power curve data: x=step sample size, y=power per test.
    const stepIndices = Array.from({ length: steps }, (_, i) => i);
    const anyKey = byKey.get(keys[0]);
    const xs = stepIndices.map(s => anyKey?.sampleSizeByStep.get(s) ?? s);
    const powerData: uPlot.AlignedData = [xs, ...keys.map(k => {
        const entry = byKey.get(k)!;
        return stepIndices.map(s => computePower(entry.pointsByStep.get(s) ?? [], ALPHA));
    })];

    let lastHoverIdx: number | null = null;
    mountPlot(powerEl, {
        width: 0, height: 0,
        title: `Power (alpha=${ALPHA})`,
        cursor: { show: true },
        scales: { x: { time: false }, y: { auto: false, range: [0, 1] } },
        axes: [
            { stroke: "#6f6354", label: "Sample size" },
            { stroke: "#6f6354", label: "Power" },
        ],
        series: [
            {},
            ...keys.map((k, i) => ({
                label: k,
                stroke: COLORS[i % COLORS.length],
                points: { show: true, size: 6 },
            })),
        ],
        hooks: {
            setCursor: [
                (u) => {
                    const idx = u.cursor.idx;
                    const next = (idx == null) ? null : idx;
                    if (next !== lastHoverIdx) {
                        lastHoverIdx = next;
                        renderSampleTable(sampleEl, groups, family, next, steps);
                    }
                },
            ],
        },
    }, powerData);

    // QQ / histogram panels: each has a small selector row and a plot area
    // that flex-fills the remaining height. Selector lets the user flip
    // between (test, step) without re-running.
    const makePanel = (host: HTMLElement, onChange: (key: string, step: number) => void) => {
        host.innerHTML = "";
        const bar = document.createElement("div");
        bar.style.cssText = "display:flex;gap:6px;margin-bottom:4px;font-size:11px;";
        const testSel = document.createElement("select");
        for (const k of keys) {
            const o = document.createElement("option");
            o.value = k; o.textContent = k;
            testSel.appendChild(o);
        }
        const stepSel = document.createElement("select");
        for (let s = 0; s < steps; s++) {
            const o = document.createElement("option");
            o.value = String(s); o.textContent = steps > 1 ? `step ${s + 1} (n=${xs[s]})` : `n=${xs[s]}`;
            stepSel.appendChild(o);
        }
        const fire = () => onChange(testSel.value, parseInt(stepSel.value, 10));
        testSel.addEventListener("change", fire);
        stepSel.addEventListener("change", fire);
        bar.appendChild(testSel);
        bar.appendChild(stepSel);
        const plot = document.createElement("div");
        plot.className = "sim-plot";
        host.appendChild(bar);
        host.appendChild(plot);
        return plot;
    };

    const qqHost = makePanel(qqEl, (k, st) => drawQQ(qqHost, byKey.get(k)!.pointsByStep.get(st) ?? [], k, st, ALPHA));
    const histHost = makePanel(histEl, (k, st) => drawHist(histHost, byKey.get(k)!.pointsByStep.get(st) ?? [], k, st));

    const initialKey = keys[0];
    drawQQ(qqHost, byKey.get(initialKey)!.pointsByStep.get(0) ?? [], initialKey, 0, ALPHA);
    drawHist(histHost, byKey.get(initialKey)!.pointsByStep.get(0) ?? [], initialKey, 0);

    // Expose for Playwright.
    (window as unknown as { __simData?: unknown }).__simData = {
        keys,
        steps,
        sampleSizes: xs,
        power: keys.map(k => {
            const e = byKey.get(k)!;
            return stepIndices.map(s => computePower(e.pointsByStep.get(s) ?? [], ALPHA));
        }),
        rowCount: rows.length,
    };
}

function drawQQ(el: HTMLElement, ps: number[], label: string, step: number, alpha: number) {
    if (ps.length === 0) return;
    const sorted = [...ps].sort((a, b) => a - b);
    const n = sorted.length;
    const xs: number[] = [], ys: number[] = [];
    for (let i = 0; i < n; i++) {
        const theoretical = (i + 1) / (n + 1);
        xs.push(-Math.log10(theoretical));
        ys.push(-Math.log10(Math.max(sorted[i], 1e-300)));
    }
    const diagY = xs.map(x => x);
    const alphaY = xs.map(() => -Math.log10(alpha));

    mountPlot(el, {
        width: 0, height: 0,
        title: `QQ — ${label} (step ${step + 1})`,
        scales: { x: { time: false }, y: { auto: true } },
        axes: [
            { stroke: "#6f6354", label: "-log10 theoretical" },
            { stroke: "#6f6354", label: "-log10 observed" },
        ],
        series: [
            {},
            { label: "obs", stroke: "#3c7ccc", points: { show: true, size: 4 }, paths: () => null },
            { label: "y=x", stroke: "#b9ab8d", dash: [4, 4] },
            { label: `α=${alpha}`, stroke: "#a62a2a", dash: [2, 4] },
        ],
    } as unknown as uPlot.Options, [xs, ys, diagY, alphaY] as unknown as uPlot.AlignedData);
}

function drawHist(el: HTMLElement, ps: number[], label: string, step: number) {
    if (ps.length === 0) return;
    const BINS = 20;
    const counts = new Array(BINS).fill(0);
    for (const p of ps) {
        const b = Math.min(BINS - 1, Math.floor(p * BINS));
        counts[b]++;
    }
    const xs = counts.map((_, i) => (i + 0.5) / BINS);
    mountPlot(el, {
        width: 0, height: 0,
        title: `Histogram — ${label} (step ${step + 1})`,
        scales: { x: { time: false, range: [0, 1] }, y: { auto: true } },
        axes: [
            { stroke: "#6f6354", label: "p-value" },
            { stroke: "#6f6354", label: "count" },
        ],
        series: [
            {},
            { label: "count", stroke: "#3c7ccc", fill: "#a8caf0",
              paths: uPlot.paths.bars!({ size: [0.9], align: 0 }) },
        ],
    }, [xs, counts] as unknown as uPlot.AlignedData);
}

// -----------------------------------------------------------------------
// HTML report export. Snapshots each plot's canvas as PNG data URI and
// stitches a self-contained, styled HTML page that opens in any browser.
// Includes summary, plots, sample-info table, and the original request
// JSON in a <details> for reproducibility.
// -----------------------------------------------------------------------
function findCanvas(host: HTMLElement): HTMLCanvasElement | null {
    return host.querySelector("canvas");
}

function escapeHtml(s: string): string {
    return s.replace(/[&<>"']/g, (c) => ({
        "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[c]!));
}

export function exportSimHtml(
    summaryEl: HTMLElement,
    powerEl: HTMLElement,
    qqEl: HTMLElement,
    histEl: HTMLElement,
    sampleEl: HTMLElement,
    request: SimRunRequest,
) {
    const power = findCanvas(powerEl);
    if (!power) throw new Error("no power plot to export");
    const qq    = findCanvas(qqEl);
    const hist  = findCanvas(histEl);

    const summary = summaryEl.textContent?.replace(/\s+/g, " ").trim() ?? "";
    const sampleHtml = sampleEl.innerHTML;
    const stamp = new Date().toISOString();

    const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>VikNGS simulation report — ${escapeHtml(stamp)}</title>
<style>
  :root { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          --bg:#f6f1e4; --panel:#fff; --text:#2d0600; --muted:#6f6354;
          --border:#d9ceb8; --accent:#3c7ccc; }
  body { margin:0; padding:24px; background:var(--bg); color:var(--text); }
  header h1 { margin:0 0 4px; font-size:20px; }
  header .stamp { color:var(--muted); font-size:12px; }
  .summary { margin:14px 0 18px; font-size:13px; color:var(--muted); }
  .summary strong { color:var(--text); }
  .panel { background:var(--panel); border:1px solid var(--border); border-radius:6px;
           padding:12px 14px; margin-bottom:14px; }
  .panel h2 { margin:0 0 10px; font-size:11px; text-transform:uppercase; letter-spacing:0.08em;
              color:var(--muted); font-weight:600; }
  img { max-width:100%; height:auto; display:block; }
  .sub-row { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
  table { width:100%; border-collapse:collapse; font-size:12px; }
  th, td { padding:4px 8px; text-align:left; border-bottom:1px solid var(--border); }
  th { color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:0.05em; font-weight:600; }
  details { margin-top:16px; font-size:12px; color:var(--muted); }
  details summary { cursor:pointer; font-size:11px; text-transform:uppercase; letter-spacing:0.06em; }
  pre { background:#fdfbf5; border:1px solid var(--border); border-radius:4px;
        padding:10px; overflow-x:auto; font-size:11px; line-height:1.4; color:var(--text); }
  footer { color:var(--muted); font-size:11px; margin-top:18px; }
</style>
</head>
<body>
<header>
  <h1>VikNGS simulation report</h1>
  <div class="stamp">Generated ${escapeHtml(stamp)}</div>
</header>
<div class="summary">${summary ? escapeHtml(summary) : ""}</div>

<section class="panel">
  <h2>Power curve</h2>
  <img src="${power.toDataURL("image/png")}" alt="power curve" />
</section>

<section class="sub-row">
  ${qq ? `<div class="panel"><h2>QQ plot</h2><img src="${qq.toDataURL("image/png")}" alt="qq" /></div>` : ""}
  ${hist ? `<div class="panel"><h2>P-value histogram</h2><img src="${hist.toDataURL("image/png")}" alt="histogram" /></div>` : ""}
</section>

<section class="panel">
  <h2>Sample info</h2>
  ${sampleHtml || "<em>(no sample info)</em>"}
</section>

<details>
  <summary>Run configuration</summary>
  <pre>${escapeHtml(JSON.stringify(request, null, 2))}</pre>
</details>

<footer>VikNGS — variant integration kit for NGS</footer>
</body>
</html>`;

    const blob = new Blob([html], { type: "text/html;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `vikngs-sim-${stamp.replace(/[:.]/g, "-")}.html`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function clearSimResults(
    summary: HTMLElement, power: HTMLElement, qq: HTMLElement, hist: HTMLElement,
    sample?: HTMLElement,
    keepSummary = false,
) {
    if (!keepSummary) summary.innerHTML = "";
    destroyAllPlots();
    for (const el of [power, qq, hist]) el.innerHTML = "";
    if (sample) sample.innerHTML = "";
}
