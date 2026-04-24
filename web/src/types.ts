// Shared message types between the UI thread and the worker.

export interface RunRequest {
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

export type WorkerMessage =
    | { kind: "log"; level: "info" | "error" | "ok"; text: string }
    | { kind: "progress"; percent: number }
    | { kind: "done"; rows: ResultRow[]; variantsParsed: number; evaluationTime: number }
    | { kind: "error"; message: string };
