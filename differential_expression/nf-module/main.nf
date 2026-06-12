#!/usr/bin/env nextflow

/*
 * bisrDE — Bulk RNA-Seq Differential Expression Pipeline (Nextflow DSL2)
 *
 * Wraps the bisrDE R package as a Nextflow process so it can sit
 * downstream of nf-core/rnaseq's merged-counts output. Phase 6.1.
 */

nextflow.enable.dsl = 2

include { BISR_DE } from './modules/local/bisr_de'

// Parameter defaults live in `nextflow.config` (the `params {}` block).
// Defining them there ensures they are visible to all process scopes, not
// just this main.nf scope.

// --------------------------------------------------------------------------
// Help
// --------------------------------------------------------------------------
def help_msg = """
bisrDE — Bulk RNA-Seq Differential Expression Pipeline (Nextflow DSL2)

Usage:
    nextflow run nf-module/main.nf \\
        --counts <FILE> \\
        --samplesheet <FILE> \\
        --runid <ID> \\
        [--outdir <DIR>] \\
        [--annotation human|mouse] \\
        [--brs_ticket BRS-XXXX] \\
        [--id_type ensembl|entrez|symbol]

Required:
    --counts        Path to merged counts.tsv (e.g. nf-core/rnaseq's rsem.merged.gene_counts.tsv).
    --samplesheet   Path to comma-delimited samplesheet.csv.
    --runid         Unique identifier for this run.

Optional:
    --outdir        Output directory.                                          [default: ${params.outdir}]
    --annotation    Genome annotation: 'mouse' or 'human'.                     [default: ${params.annotation}]
    --brs_ticket    BRS ticket identifier (e.g. BRS-1234) for the report.       [default: '']
    --id_type       Gene identifier type: 'ensembl', 'entrez', or 'symbol'.    [default: ${params.id_type}]
    --help          Show this help and exit.

Profiles (-profile <name>):
    local           Run with local R (requires bisrDE + Quarto installed locally).
    container       Run inside the dge_analysis.sif Apptainer container.
    slurm           Run on a SLURM cluster (combine with -profile container).

See nf-module/README.md for an nf-core/rnaseq integration example.
"""

if (params.help) {
    log.info help_msg
    exit 0
}

// --------------------------------------------------------------------------
// Workflow
// --------------------------------------------------------------------------
workflow {
    if (!params.counts || !params.samplesheet || !params.runid) {
        log.error("Missing required parameter(s). Run with --help for usage.")
        log.error("    counts:      ${params.counts}")
        log.error("    samplesheet: ${params.samplesheet}")
        log.error("    runid:       ${params.runid}")
        System.exit(1)
    }

    counts_ch      = Channel.fromPath(params.counts,      checkIfExists: true)
    samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)

    BISR_DE(counts_ch, samplesheet_ch)
}
