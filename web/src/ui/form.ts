// Builds the input form and returns a function that gathers field values
// into a RunRequest. Schema-driven so adding a new field is one entry.

import type { RunRequest } from "../types";

type FieldDef =
    | { id: string; label: string; type: "file"; accept?: string; optional?: boolean }
    | { id: string; label: string; type: "number"; default: number; step?: string; min?: number; max?: number }
    | { id: string; label: string; type: "text"; default: string; placeholder?: string }
    | { id: string; label: string; type: "checkbox"; default: boolean }
    | { id: string; label: string; type: "select"; default: string; options: { value: string; label: string }[] };

interface Group {
    label: string;
    fields: FieldDef[];
    collapsible?: boolean;  // render as <details>; defaults to always-visible
    startOpen?: boolean;    // only meaningful when collapsible
}

const SCHEMA: Group[] = [
    {
        label: "Files",
        fields: [
            { id: "vcf",    label: "VCF file (*.vcf)",       type: "file", accept: ".vcf" },
            { id: "sample", label: "Sample info (*.txt)",    type: "file", accept: ".txt" },
            { id: "bed",    label: "BED regions (optional)", type: "file", accept: ".bed,.txt", optional: true },
        ],
    },
    {
        label: "Test",
        collapsible: true,
        startOpen: true,
        fields: [
            { id: "statistic", label: "Statistic",        type: "select", default: "common", options: [
                { value: "common", label: "Common variant" },
                { value: "cast",   label: "Rare — CAST" },
                { value: "skat",   label: "Rare — SKAT" },
            ]},
            { id: "genotype", label: "Genotype source",   type: "select", default: "expected", options: [
                { value: "expected", label: "Expected (vRVS)" },
                { value: "call",     label: "Called (hard calls)" },
                { value: "vcf",      label: "VCF GT field" },
            ]},
            { id: "nboot",     label: "Bootstrap iterations", type: "number", default: 1, step: "1", min: 1 },
            { id: "stopEarly", label: "Stop bootstrapping early", type: "checkbox", default: false },
        ],
    },
    {
        label: "VCF filtering",
        collapsible: true,
        startOpen: false,
        fields: [
            { id: "maf",       label: "MAF threshold",            type: "number",  default: 0.05, step: "0.01", min: 0, max: 1 },
            { id: "depth",     label: "Read-depth cutoff (high/low)", type: "number", default: 30, step: "1", min: 1 },
            { id: "missing",   label: "Missing data threshold",   type: "number",  default: 0.1, step: "0.01", min: 0, max: 0.5 },
            { id: "mustPass",  label: "Require FILTER=PASS",      type: "checkbox", default: false },
            { id: "chrFilter", label: "Chromosome filter",        type: "text",    default: "", placeholder: "(all)" },
            { id: "fromPos",   label: "Min position (-1 = none)", type: "number",  default: -1, step: "1" },
            { id: "toPos",     label: "Max position (-1 = none)", type: "number",  default: -1, step: "1" },
        ],
    },
    {
        label: "Collapsing (BED only)",
        collapsible: true,
        startOpen: false,
        fields: [
            { id: "collapseMode", label: "Mode", type: "select", default: "", options: [
                { value: "",     label: "None" },
                { value: "gene", label: "By gene" },
                { value: "exon", label: "By exon" },
                { value: "k",    label: "By k variants" },
            ]},
            { id: "collapseK", label: "k (for by-k)", type: "number", default: -1, step: "1" },
        ],
    },
    {
        label: "Performance",
        collapsible: true,
        startOpen: false,
        fields: [
            { id: "batchSize", label: "Batch size", type: "number", default: 1000, step: "100", min: 1 },
        ],
    },
];

export function mountForm(root: HTMLFormElement) {
    for (const group of SCHEMA) {
        const g = group.collapsible
            ? document.createElement("details")
            : document.createElement("div");
        g.className = "group";
        if (group.collapsible) {
            if (group.startOpen) (g as HTMLDetailsElement).open = true;
            const summary = document.createElement("summary");
            summary.className = "group-label";
            summary.textContent = group.label;
            g.appendChild(summary);
        } else {
            const lbl = document.createElement("div");
            lbl.className = "group-label";
            lbl.textContent = group.label;
            g.appendChild(lbl);
        }

        for (const f of group.fields) {
            const field = document.createElement("div");
            field.className = "field";
            const label = document.createElement("label");
            label.htmlFor = `f-${f.id}`;
            label.textContent = f.label;
            field.appendChild(label);

            let input: HTMLInputElement | HTMLSelectElement;
            if (f.type === "select") {
                const sel = document.createElement("select");
                for (const opt of f.options) {
                    const o = document.createElement("option");
                    o.value = opt.value; o.textContent = opt.label;
                    if (opt.value === f.default) o.selected = true;
                    sel.appendChild(o);
                }
                input = sel;
            } else {
                const inp = document.createElement("input");
                inp.type = f.type;
                if (f.type === "checkbox") inp.checked = f.default;
                else if (f.type === "file") {
                    if (f.accept) inp.accept = f.accept;
                } else {
                    inp.value = String(f.default);
                    if (f.type === "number" && f.step) inp.step = f.step;
                    if (f.type === "number" && f.min !== undefined) inp.min = String(f.min);
                    if (f.type === "number" && f.max !== undefined) inp.max = String(f.max);
                    if (f.type === "text"   && f.placeholder)       inp.placeholder = f.placeholder;
                }
                input = inp;
            }
            input.id = `f-${f.id}`;
            input.name = f.id;
            field.appendChild(input);
            g.appendChild(field);
        }
        root.appendChild(g);
    }
}

function numField(form: HTMLFormElement, id: string): number {
    return parseFloat((form.querySelector(`#f-${id}`) as HTMLInputElement).value);
}
function textField(form: HTMLFormElement, id: string): string {
    return (form.querySelector(`#f-${id}`) as HTMLInputElement).value;
}
function boolField(form: HTMLFormElement, id: string): boolean {
    return (form.querySelector(`#f-${id}`) as HTMLInputElement).checked;
}
function selField(form: HTMLFormElement, id: string): string {
    return (form.querySelector(`#f-${id}`) as HTMLSelectElement).value;
}
function fileField(form: HTMLFormElement, id: string): File | null {
    return (form.querySelector(`#f-${id}`) as HTMLInputElement).files?.[0] ?? null;
}

export function collectForm(form: HTMLFormElement): RunRequest | { error: string } {
    const vcf = fileField(form, "vcf");
    const sample = fileField(form, "sample");
    if (!vcf)    return { error: "Please pick a VCF file." };
    if (!sample) return { error: "Please pick a sample info file." };
    return {
        vcf, sample,
        bed: fileField(form, "bed"),
        maf:       numField(form, "maf"),
        depth:     numField(form, "depth"),
        missing:   numField(form, "missing"),
        mustPass:  boolField(form, "mustPass"),
        chrFilter: textField(form, "chrFilter"),
        fromPos:   numField(form, "fromPos"),
        toPos:     numField(form, "toPos"),
        statistic: selField(form, "statistic") as "common" | "cast" | "skat",
        genotype:  selField(form, "genotype")  as "expected" | "call" | "vcf",
        nboot:     numField(form, "nboot"),
        stopEarly: boolField(form, "stopEarly"),
        collapseMode: selField(form, "collapseMode") as "" | "gene" | "exon" | "k",
        collapseK: numField(form, "collapseK"),
        batchSize: numField(form, "batchSize"),
    };
}
