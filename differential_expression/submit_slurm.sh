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
# Method knobs (v1.7): GSEA_RANK=stat|log2fc, LFC_SHRINK=apeglm|normal|none,
# INDEPENDENT_FILTERING=yes|no, e.g. to reproduce a pre-1.7 ranking:
#
#     sbatch --export=ALL,GSEA_RANK=log2fc,LFC_SHRINK=none,... submit_slurm.sh
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
GSEA_RANK="${GSEA_RANK:-stat}"                       # stat | log2fc
LFC_SHRINK="${LFC_SHRINK:-apeglm}"                   # apeglm | normal | none
INDEPENDENT_FILTERING="${INDEPENDENT_FILTERING:-no}" # yes | no (DESeq2 default is yes)

# Named in the report's Methods ("All computational analyses were performed
# on ..."). Recorded at analysis time by the pipeline; override via --export.
# (two-step default: an apostrophe inside "${VAR:-...}" breaks bash's parser)
PLATFORM_NAME="${PLATFORM_NAME:-}"
[ -n "$PLATFORM_NAME" ] || PLATFORM_NAME="VCU's High Performance Research Computing (HPRC) cluster"
export BISR_PLATFORM_NAME="$PLATFORM_NAME"

# Name of the conda environment (see environment.yml). If your site provides a
# central module instead, replace the activation block below with `module load`.
CONDA_ENV="${CONDA_ENV:-bisrde}"

# --------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------
# Locating the pipeline is not as simple as dirname "$0": sbatch COPIES this
# script into a spool directory (e.g. /var/spool/slurm/job123/slurm_script), so
# $BASH_SOURCE points there and not at the repo. Prefer $SLURM_SUBMIT_DIR (the
# directory sbatch was invoked from), fall back to the script's own location
# when run directly, and let the caller override explicitly.
if [ -n "${PIPELINE_DIR:-}" ]; then
    :                                        # caller knows best
elif [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PIPELINE_DIR="$SLURM_SUBMIT_DIR"
else
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if [ ! -f "$PIPELINE_DIR/run_analysis.sh" ]; then
    echo "ERROR: run_analysis.sh not found in: $PIPELINE_DIR" >&2
    echo "" >&2
    echo "  sbatch copies this script to a spool dir, so it locates the" >&2
    echo "  pipeline via \$SLURM_SUBMIT_DIR — submit from the" >&2
    echo "  differential_expression directory:" >&2
    echo "" >&2
    echo "    cd /path/to/bulk_rnaseq_analyses/differential_expression" >&2
    echo "    sbatch submit_slurm.sh" >&2
    echo "" >&2
    echo "  ...or pass the location explicitly:" >&2
    echo "    sbatch --export=ALL,PIPELINE_DIR=/path/to/differential_expression ..." >&2
    exit 3
fi

cd "$PIPELINE_DIR"

# Activate conda. `conda activate` needs the shell hook in a non-interactive
# job, which is why this is not just `conda activate`.
if [ -n "${CONDA_PREFIX:-}" ] && command -v Rscript >/dev/null 2>&1; then
    # Submitted with --export=ALL from an already-activated environment, so the
    # env is inherited and there is nothing to do. This is the common case, and
    # it avoids needing conda/micromamba on PATH in a non-interactive shell at
    # all (they are usually only set up by ~/.bashrc, which sbatch does not read).
    echo "using inherited environment: $CONDA_PREFIX"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
elif command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$CONDA_ENV"
else
    echo "ERROR: no activated environment, and neither conda nor micromamba is" >&2
    echo "       on PATH (sbatch does not source ~/.bashrc)." >&2
    echo "" >&2
    echo "  Easiest fix — submit from an activated environment with --export=ALL:" >&2
    echo "    micromamba activate $CONDA_ENV" >&2
    echo "    sbatch --export=ALL submit_slurm.sh" >&2
    echo "" >&2
    echo "  Or make the launcher visible to the job:" >&2
    echo "    sbatch --export=ALL,PATH=\"\$HOME/bin:\$PATH\" submit_slurm.sh" >&2
    exit 127
fi

echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-<interactive>}"
echo "env:       ${CONDA_PREFIX:-<none>}"
echo "R:         $(command -v Rscript) ($(Rscript -e 'cat(R.version.string)' 2>/dev/null))"
echo "outdir:    $OUTDIR"
echo "methods:   gsea-rank=$GSEA_RANK  lfc-shrink=$LFC_SHRINK  independent-filtering=$INDEPENDENT_FILTERING"
echo "platform:  $BISR_PLATFORM_NAME"
echo

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
args=(--counts "$COUNTS" --samplesheet "$SAMPLESHEET"
      --outdir "$OUTDIR" --runid "$RUNID"
      --annotation "$ANNOTATION" --id-type "$IDTYPE"
      --fold-change "$FOLD_CHANGE" --padj "$PADJ"
      --gsea-rank "$GSEA_RANK" --lfc-shrink "$LFC_SHRINK")
[ "$INDEPENDENT_FILTERING" = "yes" ] && args+=(--independent-filtering)
[ -n "$BRS_TICKET" ] && args+=(--brs-ticket "$BRS_TICKET")

bash run_analysis.sh "${args[@]}"

echo
echo "done. Report + logs under: $OUTDIR"
