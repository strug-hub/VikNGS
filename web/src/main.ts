import { mountForm, collectForm } from "./ui/form";
import { mountLog } from "./ui/log";
import { renderResultsTable } from "./ui/results";
import { renderManhattan } from "./ui/manhattan";
import type { RunRequest, WorkerMessage } from "./types";

const form        = document.getElementById("run-form") as HTMLFormElement;
const runBtn      = document.getElementById("run-btn") as HTMLButtonElement;
const stopBtn     = document.getElementById("stop-btn") as HTMLButtonElement;
const preloadBtn  = document.getElementById("preload-example-btn") as HTMLButtonElement;
const statusEl    = document.getElementById("run-status") as HTMLElement;
const logEl       = document.getElementById("log") as HTMLElement;
const tableEl     = document.getElementById("results-table") as HTMLElement;
const plotEl      = document.getElementById("manhattan") as HTMLElement;

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
