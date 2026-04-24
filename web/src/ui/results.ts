import type { ResultRow } from "../types";

export function renderResultsTable(el: HTMLElement, rows: ResultRow[]) {
    if (rows.length === 0) {
        el.innerHTML = "<p style='color:var(--muted);font-size:13px'>No rows passed filters.</p>";
        return;
    }
    const table = document.createElement("table");
    const thead = document.createElement("thead");
    thead.innerHTML = "<tr><th>chr</th><th>pos</th><th>ref</th><th>alt</th><th>p-value</th><th>test</th></tr>";
    table.appendChild(thead);

    const tbody = document.createElement("tbody");
    for (const r of rows) {
        const tr = document.createElement("tr");
        tr.innerHTML =
            `<td>${r.chrom}</td><td>${r.pos}</td><td>${r.ref}</td><td>${r.alt}</td>` +
            // 10-digit precision is plenty for display and round-trips cleanly
            // through parseFloat for the browser regression test.
            `<td>${r.pvalue.toPrecision(10)}</td><td>${r.testDesc}</td>`;
        tbody.appendChild(tr);
    }
    table.appendChild(tbody);
    el.innerHTML = "";
    el.appendChild(table);
}
