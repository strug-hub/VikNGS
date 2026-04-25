import { mountForm, collectForm } from "./ui/form";
import { mountSimForm, collectSimForm } from "./ui/simForm";
import { mountLog } from "./ui/log";
import { renderResultsTable } from "./ui/results";
import { renderManhattan } from "./ui/manhattan";
import { renderSimResults, clearSimResults, exportSimHtml } from "./ui/simResults";
import type { SimReportData } from "./ui/simResults";
import type { RunRequest, SimRunRequest, UiToWorker, WorkerMessage } from "./types";

const form        = document.getElementById("run-form") as HTMLFormElement;
const simForm     = document.getElementById("sim-form") as HTMLFormElement;
const runBtn      = document.getElementById("run-btn") as HTMLButtonElement;
const stopBtn     = document.getElementById("stop-btn") as HTMLButtonElement;
const simRunBtn   = document.getElementById("sim-run-btn") as HTMLButtonElement;
const simStopBtn  = document.getElementById("sim-stop-btn") as HTMLButtonElement;
const preloadBtn  = document.getElementById("preload-example-btn") as HTMLButtonElement;
const statusEl    = document.getElementById("run-status") as HTMLElement;
const simStatusEl = document.getElementById("sim-status") as HTMLElement;
const logEl       = document.getElementById("log") as HTMLElement;
const tableEl     = document.getElementById("results-table") as HTMLElement;
const plotEl      = document.getElementById("manhattan") as HTMLElement;
const simSummaryEl = document.getElementById("sim-summary") as HTMLElement;
const simPowerEl   = document.getElementById("sim-power") as HTMLElement;
const simQqEl      = document.getElementById("sim-qq") as HTMLElement;
const simHistEl    = document.getElementById("sim-hist") as HTMLElement;
const simSampleEl  = document.getElementById("sim-sample-info") as HTMLElement;
const simExportBtn = document.getElementById("sim-export-btn") as HTMLButtonElement;

mountForm(form);
mountSimForm(simForm);
const log = mountLog(logEl);

// --- Tab switching ---
function setActiveTab(tab: "analysis" | "simulation") {
    document.body.dataset.tab = tab;
    for (const b of document.querySelectorAll<HTMLElement>("nav.tabs .tab")) {
        b.classList.toggle("active", b.dataset.tab === tab);
    }
    for (const p of document.querySelectorAll<HTMLElement>(".tab-pane")) {
        p.hidden = p.dataset.tab !== tab;
    }
    // Preload example is analysis-only.
    preloadBtn.hidden = tab !== "analysis";
}

for (const b of document.querySelectorAll<HTMLButtonElement>("nav.tabs .tab")) {
    b.addEventListener("click", () => setActiveTab(b.dataset.tab as "analysis" | "simulation"));
}
setActiveTab("analysis");

// --- Worker run (shared for analysis + sim) ---
let worker: Worker | null = null;

function resetAnalysis(statusText = "") {
    runBtn.disabled = false;
    stopBtn.disabled = true;
    statusEl.className = "";
    statusEl.textContent = statusText;
}
function resetSim(statusText = "") {
    simRunBtn.disabled = false;
    simStopBtn.disabled = true;
    simStatusEl.className = "";
    simStatusEl.textContent = statusText;
}
function endRun() {
    if (worker) { worker.terminate(); worker = null; }
}

function startWorker(req: UiToWorker, onDone: (m: WorkerMessage) => void) {
    worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (ev: MessageEvent<WorkerMessage>) => {
        const m = ev.data;
        if (m.kind === "log") { log[m.level](m.text); return; }
        if (m.kind === "progress") return;
        onDone(m);
    };
    worker.onerror = (e) => {
        log.error("Worker error: " + e.message);
        resetAnalysis(); resetSim(); endRun();
    };
    worker.postMessage(req);
}

runBtn.addEventListener("click", () => {
    log.clear();
    tableEl.innerHTML = "";
    plotEl.innerHTML = "";

    const gathered = collectForm(form);
    if ("error" in gathered) { log.error(gathered.error); return; }

    runBtn.disabled = true;
    stopBtn.disabled = false;
    statusEl.className = "running";
    statusEl.textContent = "running…";
    log.info("Starting worker…");

    const req: RunRequest = gathered;
    startWorker(req, (m) => {
        if (m.kind === "done") {
            log.ok(`Done. ${m.rows.length} rows, parsed ${m.variantsParsed} variants in ${m.evaluationTime.toFixed(2)}s.`);
            renderResultsTable(tableEl, m.rows);
            renderManhattan(plotEl, m.rows);
            resetAnalysis(`done — ${m.rows.length} rows`);
            endRun();
        } else if (m.kind === "error") {
            log.error("Error: " + m.message);
            resetAnalysis("error");
            endRun();
        }
    });
});

stopBtn.addEventListener("click", () => {
    log.info("Stopping…");
    resetAnalysis("stopped");
    endRun();
});

// Last completed sim request + the data we rendered. The report exporter
// re-runs the same render in a fresh inlined HTML page so it stays
// fully interactive.
let lastSimRequest: SimRunRequest | null = null;
let lastSimReport: SimReportData | null = null;

simRunBtn.addEventListener("click", () => {
    log.clear();
    clearSimResults(simSummaryEl, simPowerEl, simQqEl, simHistEl, simSampleEl);
    simExportBtn.disabled = true;
    lastSimRequest = null;
    lastSimReport = null;

    const gathered = collectSimForm(simForm);
    if ("error" in gathered) { log.error(gathered.error); return; }

    simRunBtn.disabled = true;
    simStopBtn.disabled = false;
    simStatusEl.className = "running";
    simStatusEl.textContent = "simulating…";
    log.info("Starting simulation worker…");

    const req: SimRunRequest = gathered;
    startWorker(req, (m) => {
        if (m.kind === "sim-done") {
            log.ok(`Done. ${m.rows.length} p-values from ${m.variantsParsed} variants × ${m.steps} step(s) in ${m.evaluationTime.toFixed(2)}s.`);
            const report: SimReportData = {
                rows: m.rows,
                steps: m.steps,
                family: req.family,
                groups: req.groups,
                summary: {
                    variants: m.variantsParsed,
                    processingTime: m.processingTime,
                    evaluationTime: m.evaluationTime,
                },
            };
            renderSimResults(simSummaryEl, simPowerEl, simQqEl, simHistEl, simSampleEl, report);
            lastSimRequest = req;
            lastSimReport = report;
            simExportBtn.disabled = false;
            resetSim(`done — ${m.rows.length} p-values`);
            endRun();
        } else if (m.kind === "error") {
            log.error("Error: " + m.message);
            resetSim("error");
            endRun();
        }
    });
});

simStopBtn.addEventListener("click", () => {
    log.info("Stopping…");
    resetSim("stopped");
    endRun();
});

simExportBtn.addEventListener("click", () => {
    if (!lastSimRequest || !lastSimReport) { log.error("No simulation result to export yet."); return; }
    try {
        exportSimHtml(lastSimRequest, lastSimReport);
    } catch (e) {
        log.error("Export failed: " + (e instanceof Error ? e.message : String(e)));
    }
});

// Preload: fetch example VCF + sample info (served from /example/ via symlink
// into web/public) and populate the file inputs. Uses DataTransfer so the
// <input type="file"> actually shows the file name, not just a hidden state.
async function setFileInput(id: string, path: string, type = "text/plain") {
    const resp = await fetch(path);
    if (!resp.ok) throw new Error(`${path}: ${resp.status} ${resp.statusText}`);
    const bytes = await resp.arrayBuffer();
    const name = path.split("/").pop() ?? "file";
    const file = new File([bytes], name, { type });
    const input = document.getElementById(`f-${id}`) as HTMLInputElement;
    const dt = new DataTransfer();
    dt.items.add(file);
    input.files = dt.files;
}

preloadBtn.addEventListener("click", async () => {
    preloadBtn.disabled = true;
    try {
        log.info("Loading bundled example files…");
        await setFileInput("vcf",    "/example/example.vcf");
        await setFileInput("sample", "/example/example_info.txt");
        log.ok("Example loaded. Click Run analysis to go.");
    } catch (e) {
        log.error("Preload failed: " + (e instanceof Error ? e.message : String(e)));
    } finally {
        preloadBtn.disabled = false;
    }
});
