// Simulation form. Schema-driven like form.ts but includes an editable
// per-group table (the 5-column grid that mirrors Qt's sim_groupTbl).

import type { SimGroup, SimRunRequest } from "../types";

interface SimGroupRow extends SimGroup {
    sizeText: string;    // "200" or "200:1000" for min:max ranges
    cohort: "case" | "control";
}

const DEFAULT_ROWS: SimGroupRow[] = [
    { sizeText: "500", n: 500, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
    { sizeText: "500", n: 500, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
];

let groupRows: SimGroupRow[] = DEFAULT_ROWS.map(r => ({ ...r }));

// -----------------------------------------------------------------------
// Presets: common simulation setups that one-click-populate the form.
// -----------------------------------------------------------------------
interface Preset {
    name: string;
    title: string;
    fields: Partial<Record<
        "nsnp" | "effectSize" | "mafMin" | "mafMax" | "steps" |
        "statistic" | "collapse" | "nboot" | "stopEarly" | "seed", string | number | boolean
    >>;
    groups: SimGroupRow[];
}

const PRESETS: Preset[] = [
    {
        name: "Null (Type I)",
        title: "Common-variant null — 500 vs 500, effect=1.0",
        fields: { nsnp: 500, effectSize: 1.0, mafMin: 0.05, mafMax: 0.5, steps: 1, statistic: "common", collapse: 1, nboot: 1, stopEarly: false, seed: 0 },
        groups: [
            { sizeText: "500", n: 500, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "500", n: 500, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
        ],
    },
    {
        name: "Common effect",
        title: "Common-variant OR=1.4, 400 vs 400, high depth",
        fields: { nsnp: 300, effectSize: 1.4, mafMin: 0.05, mafMax: 0.5, steps: 1, statistic: "common", collapse: 1, nboot: 1, stopEarly: false, seed: 0 },
        groups: [
            { sizeText: "400", n: 400, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "400", n: 400, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
        ],
    },
    {
        name: "Power sweep",
        title: "Common-variant OR=1.3, sample size 200→2000 in 5 steps",
        fields: { nsnp: 200, effectSize: 1.3, mafMin: 0.05, mafMax: 0.5, steps: 5, statistic: "common", collapse: 1, nboot: 1, stopEarly: false, seed: 0 },
        groups: [
            { sizeText: "200:2000", n: 200, nIncrement: 450, isCase: true,  cohort: "case",    meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "200:2000", n: 200, nIncrement: 450, isCase: false, cohort: "control", meanDepth: 20.0, sdDepth: 2.0, errorRate: 0.01, readDepth: "high" },
        ],
    },
    {
        name: "Rare SKAT",
        title: "Rare-variant SKAT, 500 variants in 25-variant gene sets, OR=2.5",
        fields: { nsnp: 500, effectSize: 2.5, mafMin: 0.001, mafMax: 0.05, steps: 1, statistic: "skat", collapse: 25, nboot: 1000, stopEarly: true, seed: 0 },
        groups: [
            { sizeText: "1000", n: 1000, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "1000", n: 1000, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
        ],
    },
    {
        name: "Rare CAST",
        title: "Rare-variant CAST, 500 variants × 25-variant sets, OR=2.5",
        fields: { nsnp: 500, effectSize: 2.5, mafMin: 0.001, mafMax: 0.05, steps: 1, statistic: "cast", collapse: 25, nboot: 1000, stopEarly: true, seed: 0 },
        groups: [
            { sizeText: "1000", n: 1000, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "1000", n: 1000, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
        ],
    },
    {
        name: "Depth imbalance",
        title: "Unequal depth — cases high (30×) vs controls low (5×). Tests vRVS.",
        fields: { nsnp: 300, effectSize: 1.0, mafMin: 0.05, mafMax: 0.5, steps: 1, statistic: "common", collapse: 1, nboot: 1, stopEarly: false, seed: 0 },
        groups: [
            { sizeText: "500", n: 500, nIncrement: 0, isCase: true,  cohort: "case",    meanDepth: 30.0, sdDepth: 3.0, errorRate: 0.01, readDepth: "high" },
            { sizeText: "500", n: 500, nIncrement: 0, isCase: false, cohort: "control", meanDepth: 5.0,  sdDepth: 1.5, errorRate: 0.01, readDepth: "low"  },
        ],
    },
];

let rerenderGroupTableRef: (() => void) | null = null;

function applyPreset(root: HTMLFormElement, preset: Preset) {
    for (const [id, val] of Object.entries(preset.fields)) {
        const inp = root.querySelector<HTMLInputElement | HTMLSelectElement>(`#sf-${id}`);
        if (!inp) continue;
        if (inp.type === "checkbox") {
            (inp as HTMLInputElement).checked = !!val;
        } else {
            inp.value = String(val);
        }
        // Fire change so any listeners update.
        inp.dispatchEvent(new Event("change", { bubbles: true }));
    }
    groupRows = preset.groups.map(r => ({ ...r }));
    rerenderGroupTableRef?.();
}

function parseSize(text: string): { n: number; nIncrement: number } {
    const m = text.match(/^(\d+)\s*:\s*(\d+)$/);
    if (m) {
        const lo = parseInt(m[1], 10);
        const hi = parseInt(m[2], 10);
        return { n: lo, nIncrement: hi - lo };
    }
    return { n: parseInt(text, 10) || 0, nIncrement: 0 };
}

function rerenderGroupTable(tbody: HTMLTableSectionElement) {
    tbody.innerHTML = "";
    groupRows.forEach((row, idx) => {
        const tr = document.createElement("tr");

        const tdSize = document.createElement("td");
        const sizeInp = document.createElement("input");
        sizeInp.type = "text";
        sizeInp.value = row.sizeText;
        sizeInp.placeholder = "e.g. 500 or 500:2000";
        sizeInp.addEventListener("input", () => {
            row.sizeText = sizeInp.value;
            const { n, nIncrement } = parseSize(sizeInp.value);
            row.n = n; row.nIncrement = nIncrement;
        });
        tdSize.appendChild(sizeInp);

        const tdCohort = document.createElement("td");
        const cohortSel = document.createElement("select");
        for (const opt of ["case", "control"]) {
            const o = document.createElement("option");
            o.value = opt; o.textContent = opt;
            if (opt === row.cohort) o.selected = true;
            cohortSel.appendChild(o);
        }
        cohortSel.addEventListener("change", () => {
            row.cohort = cohortSel.value as "case" | "control";
            row.isCase = row.cohort === "case";
        });
        tdCohort.appendChild(cohortSel);

        const makeNumTd = (get: () => number, set: (v: number) => void, step = "0.1") => {
            const td = document.createElement("td");
            const inp = document.createElement("input");
            inp.type = "number"; inp.step = step; inp.value = String(get());
            inp.addEventListener("input", () => set(parseFloat(inp.value)));
            td.appendChild(inp);
            return td;
        };

        const tdDepth = makeNumTd(() => row.meanDepth, v => row.meanDepth = v, "1");
        const tdSd    = makeNumTd(() => row.sdDepth,   v => row.sdDepth = v,   "0.1");
        const tdErr   = makeNumTd(() => row.errorRate, v => row.errorRate = v, "0.001");

        const tdHigh = document.createElement("td");
        const highSel = document.createElement("select");
        for (const opt of ["high", "low"]) {
            const o = document.createElement("option");
            o.value = opt; o.textContent = opt;
            if (opt === row.readDepth) o.selected = true;
            highSel.appendChild(o);
        }
        highSel.addEventListener("change", () => { row.readDepth = highSel.value as "high" | "low"; });
        tdHigh.appendChild(highSel);

        const tdDel = document.createElement("td");
        tdDel.className = "row-del";
        const delBtn = document.createElement("button");
        delBtn.type = "button"; delBtn.textContent = "×"; delBtn.className = "secondary";
        delBtn.title = "Remove this group";
        delBtn.addEventListener("click", () => {
            groupRows.splice(idx, 1);
            rerenderGroupTable(tbody);
        });
        tdDel.appendChild(delBtn);

        tr.appendChild(tdSize);
        tr.appendChild(tdCohort);
        tr.appendChild(tdDepth);
        tr.appendChild(tdSd);
        tr.appendChild(tdErr);
        tr.appendChild(tdHigh);
        tr.appendChild(tdDel);
        tbody.appendChild(tr);
    });
}

export function mountSimForm(root: HTMLFormElement) {
    root.innerHTML = "";

    // Preset strip across the top of the form.
    const presets = document.createElement("div");
    presets.className = "sim-presets";
    const plabel = document.createElement("span");
    plabel.className = "presets-label";
    plabel.textContent = "Presets:";
    presets.appendChild(plabel);
    for (const p of PRESETS) {
        const b = document.createElement("button");
        b.type = "button";
        b.textContent = p.name;
        b.title = p.title;
        b.addEventListener("click", () => applyPreset(root, p));
        presets.appendChild(b);
    }
    root.appendChild(presets);

    const mkGroup = (label: string, open = true) => {
        const g = document.createElement("details");
        g.className = "group";
        if (open) g.open = true;
        const s = document.createElement("summary");
        s.className = "group-label";
        s.textContent = label;
        g.appendChild(s);
        return g;
    };

    const mkNumField = (id: string, label: string, def: number, step = "any", min?: number, max?: number) => {
        const field = document.createElement("div");
        field.className = "field";
        const lbl = document.createElement("label");
        lbl.htmlFor = `sf-${id}`; lbl.textContent = label;
        field.appendChild(lbl);
        const inp = document.createElement("input");
        inp.type = "number"; inp.id = `sf-${id}`; inp.name = id;
        inp.step = step; inp.value = String(def);
        if (min !== undefined) inp.min = String(min);
        if (max !== undefined) inp.max = String(max);
        field.appendChild(inp);
        return field;
    };

    const mkSelectField = (id: string, label: string, def: string, opts: { value: string; label: string }[]) => {
        const field = document.createElement("div");
        field.className = "field";
        const lbl = document.createElement("label");
        lbl.htmlFor = `sf-${id}`; lbl.textContent = label;
        field.appendChild(lbl);
        const sel = document.createElement("select");
        sel.id = `sf-${id}`; sel.name = id;
        for (const o of opts) {
            const oo = document.createElement("option");
            oo.value = o.value; oo.textContent = o.label;
            if (o.value === def) oo.selected = true;
            sel.appendChild(oo);
        }
        field.appendChild(sel);
        return field;
    };

    const mkCheckField = (id: string, label: string, def: boolean) => {
        const field = document.createElement("div");
        field.className = "field";
        const lbl = document.createElement("label");
        lbl.htmlFor = `sf-${id}`; lbl.textContent = label;
        field.appendChild(lbl);
        const inp = document.createElement("input");
        inp.type = "checkbox"; inp.id = `sf-${id}`; inp.name = id; inp.checked = def;
        field.appendChild(inp);
        return field;
    };

    // Variants group
    const gVar = mkGroup("Variants", true);
    gVar.appendChild(mkNumField("nsnp", "Number of variants", 200, "1", 1));
    gVar.appendChild(mkNumField("effectSize", "Odds ratio (1 = null)", 1.0, "0.05", 0));
    gVar.appendChild(mkNumField("mafMin", "MAF min", 0.05, "0.01", 0.0001, 0.5));
    gVar.appendChild(mkNumField("mafMax", "MAF max", 0.5, "0.01", 0.0001, 0.5));
    root.appendChild(gVar);

    // Groups table
    const gGroups = mkGroup("Groups (case/control)", true);
    const table = document.createElement("table");
    table.className = "sim-groups";
    const thead = document.createElement("thead");
    thead.innerHTML = `<tr>
        <th>Sample size<br><span style="font-weight:400;font-size:10px;">n or min:max</span></th>
        <th>Cohort</th>
        <th>Mean depth</th>
        <th>Depth SD</th>
        <th>Error</th>
        <th>Depth cat.</th>
        <th></th>
    </tr>`;
    table.appendChild(thead);
    const tbody = document.createElement("tbody");
    table.appendChild(tbody);
    gGroups.appendChild(table);
    const rowActions = document.createElement("div");
    rowActions.className = "row-actions";
    const addBtn = document.createElement("button");
    addBtn.type = "button"; addBtn.textContent = "+ Add group"; addBtn.className = "secondary";
    addBtn.addEventListener("click", () => {
        groupRows.push({
            sizeText: "500", n: 500, nIncrement: 0,
            isCase: false, cohort: "control",
            meanDepth: 20, sdDepth: 2, errorRate: 0.01, readDepth: "high",
        });
        rerenderGroupTable(tbody);
    });
    rowActions.appendChild(addBtn);
    gGroups.appendChild(rowActions);
    rerenderGroupTable(tbody);
    rerenderGroupTableRef = () => rerenderGroupTable(tbody);
    root.appendChild(gGroups);

    // Power steps
    const gSteps = mkGroup("Power steps", false);
    gSteps.appendChild(mkNumField("steps", "Steps (1 = no power sweep)", 1, "1", 1));
    root.appendChild(gSteps);

    // Test
    const gTest = mkGroup("Test", false);
    gTest.appendChild(mkSelectField("statistic", "Statistic", "common", [
        { value: "common", label: "Common variant" },
        { value: "cast",   label: "Rare — CAST" },
        { value: "skat",   label: "Rare — SKAT" },
        { value: "calpha", label: "Rare — C-alpha" },
    ]));
    gTest.appendChild(mkNumField("collapse", "Collapse size (rare)", 1, "1", 1));
    gTest.appendChild(mkNumField("nboot", "Bootstrap iterations (1 = asymptotic)", 1, "1", 1));
    gTest.appendChild(mkCheckField("stopEarly", "Stop bootstrapping early", false));
    root.appendChild(gTest);

    // Advanced
    const gAdv = mkGroup("Advanced", false);
    gAdv.appendChild(mkNumField("seed", "RNG seed (0 = nondeterministic)", 0, "1", 0));
    root.appendChild(gAdv);
}

function n(form: HTMLFormElement, id: string): number {
    return parseFloat((form.querySelector(`#sf-${id}`) as HTMLInputElement).value);
}
function b(form: HTMLFormElement, id: string): boolean {
    return (form.querySelector(`#sf-${id}`) as HTMLInputElement).checked;
}
function s(form: HTMLFormElement, id: string): string {
    return (form.querySelector(`#sf-${id}`) as HTMLSelectElement).value;
}

export function collectSimForm(form: HTMLFormElement): SimRunRequest | { error: string } {
    if (groupRows.length < 2) return { error: "Binomial simulation needs at least one case and one control group." };
    const cases = groupRows.filter(r => r.isCase).length;
    const ctrls = groupRows.filter(r => !r.isCase).length;
    if (cases < 1 || ctrls < 1) return { error: "Need at least one case and one control group." };

    const steps = Math.max(1, Math.floor(n(form, "steps")));
    const seed = Math.floor(n(form, "seed"));

    const groups: SimGroup[] = groupRows.map(r => ({
        n: r.n,
        nIncrement: r.nIncrement,
        isCase: r.isCase,
        meanDepth: r.meanDepth,
        sdDepth: r.sdDepth,
        errorRate: r.errorRate,
        readDepth: r.readDepth,
    }));

    return {
        kind: "sim",
        nsnp: Math.floor(n(form, "nsnp")),
        effectSize: n(form, "effectSize"),
        mafMin: n(form, "mafMin"),
        mafMax: n(form, "mafMax"),
        steps,
        family: "binomial",
        statistic: s(form, "statistic") as "common" | "cast" | "skat" | "calpha",
        collapse: Math.floor(n(form, "collapse")),
        nboot: Math.floor(n(form, "nboot")),
        stopEarly: b(form, "stopEarly"),
        groups,
        seed,
    };
}
