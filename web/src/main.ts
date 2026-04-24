import { mountForm, collectForm } from "./ui/form";
import { mountLog } from "./ui/log";
import { renderResultsTable } from "./ui/results";
import { renderManhattan } from "./ui/manhattan";
import type { RunRequest, WorkerMessage } from "./types";

const form      = document.getElementById("run-form") as HTMLFormElement;
const runBtn    = document.getElementById("run-btn") as HTMLButtonElement;
const stopBtn   = document.getElementById("stop-btn") as HTMLButtonElement;
const statusEl  = document.getElementById("run-status") as HTMLElement;
const logEl     = document.getElementById("log") as HTMLElement;
const tableEl   = document.getElementById("results-table") as HTMLElement;
const plotEl    = document.getElementById("manhattan") as HTMLElement;

mountForm(form);
const log = mountLog(logEl);

let worker: Worker | null = null;

function endRun(statusText = "") {
    runBtn.disabled = false;
    stopBtn.disabled = true;
    statusEl.className = "";
    statusEl.textContent = statusText;
    if (worker) { worker.terminate(); worker = null; }
}

runBtn.addEventListener("click", async () => {
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

    worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (ev: MessageEvent<WorkerMessage>) => {
        const m = ev.data;
        switch (m.kind) {
            case "log":      log[m.level](m.text); break;
            case "progress": break;  // Phase G
            case "done":
                log.ok(`Done. ${m.rows.length} rows, parsed ${m.variantsParsed} variants in ${m.evaluationTime.toFixed(2)}s.`);
                renderResultsTable(tableEl, m.rows);
                renderManhattan(plotEl, m.rows);
                endRun(`done — ${m.rows.length} rows`);
                break;
            case "error":
                log.error("Error: " + m.message);
                endRun("error");
                break;
        }
    };
    worker.onerror = (e) => {
        log.error("Worker error: " + e.message);
        endRun();
    };

    const req: RunRequest = gathered;
    worker.postMessage(req);
});

stopBtn.addEventListener("click", () => {
    log.info("Stopping…");
    endRun("stopped");
});
