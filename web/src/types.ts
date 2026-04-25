// Shared message types between the UI thread and the worker.

export interface RunRequest {
    kind?: "analysis";
    vcf: File;
    sample: File;
    bed: File | null;
    maf: number;
    depth: number;
    missing: number;
    mustPass: boolean;
    chrFilter: string;
    fromPos: number;
    toPos: number;
    statistic: "common" | "cast" | "skat";
    genotype: "expected" | "call" | "vcf";
    nboot: number;
    stopEarly: boolean;
    collapseMode: "" | "gene" | "exon" | "k";
    collapseK: number;
    batchSize: number;
}

export interface ResultRow {
    chrom: string;
    pos: number;
    ref: string;
    alt: string;
    pvalue: number;
    testDesc: string;
}

export interface SimGroup {
    n: number;
    nIncrement: number;
    isCase: boolean;          // ignored when family=normal
    normalMean: number;       // only used when family=normal
    normalSd: number;         // only used when family=normal
    meanDepth: number;
    sdDepth: number;
    errorRate: number;
    readDepth: "high" | "low";
}

export interface SimRunRequest {
    kind: "sim";
    nsnp: number;
    effectSize: number;
    mafMin: number;
    mafMax: number;
    steps: number;
    family: "binomial" | "normal";
    statistic: "common" | "cast" | "skat" | "calpha";
    collapse: number;
    nboot: number;
    stopEarly: boolean;
    covariate: number;        // <0 disables; otherwise correlation strength
    corX: boolean;            // true = correlate with X (genotype)
    groups: SimGroup[];
    seed: number;
}

export interface SimResultRow {
    stepIdx: number;
    sampleSize: number;
    testIdx: number;
    statName: string;
    genotypeSource: string;  // "true" | "expected" | "call"
    variantIdx: number;
    pvalue: number;
}

export interface SimDone {
    kind: "sim-done";
    rows: SimResultRow[];
    steps: number;
    processingTime: number;
    evaluationTime: number;
    variantsParsed: number;
}

export type WorkerMessage =
    | { kind: "log"; level: "info" | "error" | "ok"; text: string }
    | { kind: "progress"; percent: number }
    | { kind: "done"; rows: ResultRow[]; variantsParsed: number; evaluationTime: number }
    | SimDone
    | { kind: "error"; message: string };

// Discriminated union for messages UI → worker. Analysis is the legacy
// request (no `kind` on the wire for backwards-compat; worker type-narrows).
export type UiToWorker = RunRequest | SimRunRequest;
