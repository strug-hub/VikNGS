// Per-sample drill-down panel. Renders the table of sample rows for one
// variant (called when the user clicks a Manhattan point).

import type { AnalysisDetail, SampleGenotype } from "../types";

function fmtDosage(v: number): string {
    return Number.isNaN(v) ? "—" : v.toFixed(3);
}

function cohort(s: SampleGenotype): "case" | "control" | "" {
    if (Number.isNaN(s.phenotype)) return "";
    return s.phenotype > 0.5 ? "case" : "control";
}

export function renderDetail(host: HTMLElement, detail: AnalysisDetail) {
    if (detail.errorMessage) {
        host.innerHTML = `<div style="color:var(--error)">${detail.errorMessage}</div>`;
        return;
    }

    // Decide which dosage columns to show — only those with at least one
    // non-NaN value in the result.
    const samples = detail.samples;
    const hasTrue = samples.some(s => !Number.isNaN(s.trueDosage));
    const hasExp  = samples.some(s => !Number.isNaN(s.expectedDosage));
    const hasCall = samples.some(s => !Number.isNaN(s.callDosage));
    const hasVcf  = samples.some(s => !Number.isNaN(s.vcfDosage));

    let html = '<table><thead><tr><th>#</th><th>cohort</th><th>group</th><th>y</th>';
    if (hasTrue) html += '<th>true GT</th>';
    if (hasExp)  html += '<th>expected</th>';
    if (hasCall) html += '<th>call</th>';
    if (hasVcf)  html += '<th>VCF GT</th>';
    html += '</tr></thead><tbody>';

    for (const s of samples) {
        const c = cohort(s);
        const phen = Number.isNaN(s.phenotype) ? "—" : s.phenotype.toFixed(3);
        html += `<tr><td>${s.sampleIdx}</td><td class="${c}">${c || "—"}</td><td>${s.group}</td><td>${phen}</td>`;
        if (hasTrue) html += `<td class="${Number.isNaN(s.trueDosage) ? "missing" : ""}">${fmtDosage(s.trueDosage)}</td>`;
        if (hasExp)  html += `<td class="${Number.isNaN(s.expectedDosage) ? "missing" : ""}">${fmtDosage(s.expectedDosage)}</td>`;
        if (hasCall) html += `<td class="${Number.isNaN(s.callDosage) ? "missing" : ""}">${fmtDosage(s.callDosage)}</td>`;
        if (hasVcf)  html += `<td class="${Number.isNaN(s.vcfDosage) ? "missing" : ""}">${fmtDosage(s.vcfDosage)}</td>`;
        html += '</tr>';
    }
    html += '</tbody></table>';
    host.innerHTML = html;
}
