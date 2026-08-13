#!/usr/bin/env bash
#
#SBATCH --job-name=bisrde
#SBATCH --output=slurm-bisrde-%j.out
#SBATCH --error=slurm-bisrde-%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#
# submit_slurm.sh — run the bisrDE pipeline as a batch job.
#
# Sizing: DESeq2 + the four GSEA backends are memory-hungry rather than
# CPU-hungry; 32G/4cpu comfortably handles a typical 20-40 sample human or
# mouse run. Raise --mem before --cpus-per-task if a run is killed.
#
# Usage — edit the CONFIG block, then:
#
#     sbatch submit_slurm.sh
#
# or override without editing:
#
#     sbatch --export=ALL,COUNTS=/path/counts.tsv,SAMPLESHEET=/path/ss.csv,\
#     OUTDIR=/path/results,RUNID=my_run,ANNOTATION=human submit_slurm.sh
#
# Check on it with:  squeue -u $USER      /  tail -f slurm-bisrde-<jobid>.out

set -euo pipefail

# --------------------------------------------------------------------------
# CONFIG — edit these, or pass them via --export (values below are fallbacks).
# --------------------------------------------------------------------------
COUNTS="${COUNTS:-/path/to/salmon.merged.gene_counts.tsv}"
SAMPLESHEET="${SAMPLESHEET:-/path/to/samplesheet.csv}"
OUTDIR="${OUTDIR:-$PWD/results}"
RUNID="${RUNID:-run_$(date +%Y%m%d_%H%M%S)}"
ANNOTATION="${ANNOTATION:-human}"          # human | mouse
IDTYPE="${IDTYPE:-ensembl}"                # ensembl | entrez | symbol
BRS_TICKET="${BRS_TICKET:-}"               # e.g. BRS-1234 (optional)
FOLD_CHANGE="${FOLD_CHANGE:-1.5}"
PADJ="${PADJ:-0.05}"

# Name of the conda environment (see environment.yml). If your site provides a
# central module instead, replace the activation block below with `module load`.
CONDA_ENV="${CONDA_ENV:-bisrde}"

# --------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

# Activate conda. `conda activate` needs the shell hook in a non-interactive
# job, which is why this is not just `conda activate`.
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
elif command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$CONDA_ENV"
else
    echo "ERROR: neither conda nor micromamba found on PATH." >&2
    echo "See the HPC deployment section of README.md." >&2
    exit 127
fi

echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-<interactive>}"
echo "env:       ${CONDA_PREFIX:-<none>}"
echo "R:         $(command -v Rscript) ($(Rscript -e 'cat(R.version.string)' 2>/dev/null))"
echo "outdir:    $OUTDIR"
echo

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
args=(--counts "$COUNTS" --samplesheet "$SAMPLESHEET"
      --outdir "$OUTDIR" --runid "$RUNID"
      --annotation "$ANNOTATION" --id-type "$IDTYPE"
      --fold-change "$FOLD_CHANGE" --padj "$PADJ")
[ -n "$BRS_TICKET" ] && args+=(--brs-ticket "$BRS_TICKET")

bash run_analysis.sh "${args[@]}"

echo
echo "done. Report + logs under: $OUTDIR"
