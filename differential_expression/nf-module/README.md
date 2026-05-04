# bisrDE Nextflow module

DSL2 wrapper around the [`bisrDE` R package](../bisrDE/) for VCU Massey BISR's
bulk RNA-seq differential expression pipeline. Designed to drop in downstream
of [`nf-core/rnaseq`](https://nf-co.re/rnaseq/)'s merged-counts output.

## Requirements

- Nextflow >= 22.10 (uses DSL2). Check via `nextflow -version`.
- Java 11–21. Note: Nextflow 23.x has been known to fail on Java 25; if you
  see "Cannot find Java", point `JAVA_HOME` at a 21-or-earlier JDK:
  ```sh
  export JAVA_HOME=$(/usr/libexec/java_home -v 21)
  ```
- For the `local` profile: R 4.2+ with the `bisrDE` package available
  (installed, or accessible via `devtools::load_all()` from `../bisrDE/`),
  plus a [Quarto CLI](https://quarto.org) install for `bisrDE::generate_report`.
- For the `container` profile: `apptainer` (or `singularity`) plus the
  `dge_analysis.sif` image built from `../dge_analysis.def` (Phase 8).

## Quickstart

Run on the bundled mouse example fixture (10 genes × 6 samples):

```sh
cd /path/to/differential_expression
nextflow run nf-module/main.nf \
    --counts        assets/example_counts.tsv \
    --samplesheet   assets/example_samplesheet.csv \
    --runid         demo_run \
    --annotation    mouse \
    --outdir        results/
```

Output:
- `results/rnaseq_analysis_<timestamp>.html` — the main deliverable
- `results/analysis_results_<timestamp>.rds` — re-renderable analysis bundle
- `results/data/`, `results/figures/`, `results/logs/` — supporting artefacts
- `results/nextflow_report.html`, `results/nextflow_trace.txt` — Nextflow run report

Get help:

```sh
nextflow run nf-module/main.nf --help
```

## Parameters

| Param           | Required | Default     | Description                                          |
|-----------------|----------|-------------|------------------------------------------------------|
| `--counts`      | yes      | —           | Merged counts TSV (e.g. `rsem.merged.gene_counts.tsv`) |
| `--samplesheet` | yes      | —           | Comma-delimited samplesheet (`SampleID, GroupID, ...contrasts`) |
| `--runid`       | yes      | —           | Unique identifier for this run                       |
| `--outdir`      | no       | `./results` | Output directory                                     |
| `--annotation`  | no       | `mouse`     | `mouse` or `human`                                   |
| `--brs_ticket`  | no       | `''`        | BRS ticket identifier (rendered as report subtitle)  |
| `--id_type`     | no       | `ensembl`   | `ensembl`, `entrez`, or `symbol`                     |

## Profiles

```sh
nextflow run nf-module/main.nf -profile <name> ...
```

- `local` — run with the host R environment. Requires `bisrDE` + Quarto installed locally.
- `container` — run inside `dge_analysis.sif` (auto-bound via Apptainer). Phase 8 rebuilds this container.
- `slurm` — submit as a SLURM job. Combine with `container`: `-profile container,slurm`.

## Integration with `nf-core/rnaseq`

The bisrDE module consumes the `merged_gene_counts` channel emitted by
`nf-core/rnaseq` and a samplesheet that describes the contrasts. Wire it
into a parent workflow like this:

```groovy
// parent.nf
nextflow.enable.dsl = 2

include { RNASEQ }  from 'nf-core/rnaseq/main'
include { BISR_DE } from './nf-module/modules/local/bisr_de'

workflow {
    // 1. Run nf-core/rnaseq up through merged-counts
    RNASEQ()

    // 2. Take the merged gene-counts TSV from nf-core/rnaseq's output
    //    (named output channel `merged_gene_counts` on the salmon/star paths;
    //    the path may be slightly different per nf-core/rnaseq version —
    //    check the pipeline's emit declarations).
    counts_ch = RNASEQ.out.merged_gene_counts

    // 3. The samplesheet for DE is your design matrix, NOT nf-core/rnaseq's
    //    input samplesheet. Build it separately with SampleID, GroupID, and
    //    one column per contrast (1 = exp, 0 = ctrl, blank = exclude).
    samplesheet_ch = Channel.fromPath(params.de_samplesheet, checkIfExists: true)

    // 4. Run bisrDE
    BISR_DE(counts_ch, samplesheet_ch)
}
```

Run with:

```sh
nextflow run parent.nf \
    --input            <rnaseq_samplesheet.csv> \
    --de_samplesheet   <de_samplesheet.csv> \
    --runid            project_X \
    --annotation       human \
    --brs_ticket       BRS-1234 \
    --outdir           results/ \
    -profile container
```

## Testing

Unit + smoke tests are written in [`nf-test`](https://www.nf-test.com).
Install via `curl -fsSL https://code.askimed.com/install/nf-test | bash` then:

```sh
cd nf-module/
nf-test test tests/main.nf.test
```

The test runs the BISR_DE process against `assets/example_*` and asserts
that `report`, `rds`, `figures/`, and `data/` outputs are produced.
Wall-clock: ~90 s on a M-series Mac.
