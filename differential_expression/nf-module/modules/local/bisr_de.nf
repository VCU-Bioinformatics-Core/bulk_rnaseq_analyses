/*
 * BISR_DE — single per-run process that delegates to the parent's de.R
 * thin wrapper, which in turn calls bisrDE::run_pipeline +
 * bisrDE::generate_report. Phase 6.1.
 *
 * Inputs:
 *   - counts: merged gene-counts TSV (e.g. nf-core/rnaseq's
 *     rsem.merged.gene_counts.tsv).
 *   - samplesheet: comma-delimited samplesheet.csv with SampleID, GroupID,
 *     and one column per contrast.
 *
 * Outputs (published to params.outdir via publishDir):
 *   - report:  rnaseq_analysis_*.html (the main deliverable).
 *   - rds:     analysis_results_*.rds (analysis bundle for re-rendering).
 *   - data:    data/ subdir with DESeq2 + GSEA CSVs.
 *   - figures: figures/ subdir with volcano, heatmap, PCA, QC PNGs.
 *   - logs:    logs/ subdir with the per-run session log.
 */

process BISR_DE {
    tag { "${params.runid}" }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path counts
    path samplesheet

    output:
    path "rnaseq_analysis_*.html",  emit: report
    path "analysis_results_*.rds",  emit: rds
    path "data",                     emit: data
    path "figures",                  emit: figures
    path "logs",                     emit: logs

    script:
    def brs_arg = params.brs_ticket?.trim() ? "--brs-ticket ${params.brs_ticket}" : ""
    // v1.5.0 — optional sample/group/contrast selection flags.
    def excl_samp_arg = params.exclude_samples?.trim()   ? "--exclude-samples ${params.exclude_samples}"     : ""
    def excl_grp_arg  = params.exclude_groups?.trim()    ? "--exclude-groups ${params.exclude_groups}"       : ""
    def incl_con_arg  = params.include_contrasts?.trim() ? "--include-contrasts ${params.include_contrasts}" : ""
    def excl_con_arg  = params.exclude_contrasts?.trim() ? "--exclude-contrasts ${params.exclude_contrasts}" : ""
    def de_dir  = "${projectDir}/.."
    """
    # Capture absolute paths BEFORE cd-ing away (Nextflow stages inputs as
    # symlinks in the work dir).
    COUNTS_ABS=\$(readlink -f ${counts} 2>/dev/null || realpath ${counts})
    SS_ABS=\$(readlink -f ${samplesheet} 2>/dev/null || realpath ${samplesheet})
    OUT_ABS=\$PWD

    # cd into the parent project dir so R picks up its .Rprofile (which
    # activates renv and exposes the bisrDE package via devtools::load_all).
    cd ${de_dir}

    Rscript de.R \\
        --counts \$COUNTS_ABS \\
        --samplesheet \$SS_ABS \\
        --outdir \$OUT_ABS \\
        --runid ${params.runid} \\
        --annotation ${params.annotation} \\
        --id-type ${params.id_type} \\
        ${brs_arg} ${excl_samp_arg} ${excl_grp_arg} ${incl_con_arg} ${excl_con_arg}
    """
}
