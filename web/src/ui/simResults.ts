// Simulation result rendering: power curve, QQ plot, histogram.
// Mirrors SimPlotWindowPlotter.cpp (Qt). Plots fill available vertical
// space and re-fit on container resize.

import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";
import type { SimResultRow } from "../types";

const COLORS = ["#3c7ccc", "#a62a2a", "#3f7a4a", "#b68a00", "#8c46c6", "#2d8ba3"];

export interface SimRenderInput {
    rows: SimResultRow[];
    steps: number;
    summary: { variants: number; processingTime: number; evaluationTime: number };
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
    input: SimRenderInput,
) {
    const { rows, steps, summary } = input;
    summaryEl.innerHTML = `
      <strong>${summary.variants}</strong> variants &middot;
      processing <strong>${summary.processingTime.toFixed(2)}s</strong> &middot;
      evaluation <strong>${summary.evaluationTime.toFixed(2)}s</strong> &middot;
      <strong>${rows.length}</strong> p-values
    `;
    clearSimResults(summaryEl, powerEl, qqEl, histEl, /*keepSummary=*/true);

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

export function clearSimResults(
    summary: HTMLElement, power: HTMLElement, qq: HTMLElement, hist: HTMLElement,
    keepSummary = false,
) {
    if (!keepSummary) summary.innerHTML = "";
    destroyAllPlots();
    for (const el of [power, qq, hist]) el.innerHTML = "";
}
