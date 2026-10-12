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
# Method knobs (defaults: GSEA_RANK=stat, LFC_SHRINK=apeglm, INDEPENDENT_FILTERING=yes),
# e.g. to reproduce a pre-1.7 ranking and the pre-1.9 padj behaviour:
#
#     sbatch --export=ALL,GSEA_RANK=log2fc,LFC_SHRINK=none,INDEPENDENT_FILTERING=no,... submit_slurm.sh
#
# Runtime (default RUNTIME=auto: the container when apptainer is on PATH and the
# image has been pulled, otherwise conda). To use the prebuilt image, pull it
# once on a login node with `bash pull_container.sh`, then:
#
#     sbatch --export=ALL,RUNTIME=container,BISR_SIF_DIR=/shared/bisrde,... submit_slurm.sh
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
ANALYST="${ANALYST:-}"                     # name shown in the report header (optional, default Mikail Bala)
FOLD_CHANGE="${FOLD_CHANGE:-1.5}"
PADJ="${PADJ:-0.05}"
GSEA_RANK="${GSEA_RANK:-stat}"                       # stat | log2fc
LFC_SHRINK="${LFC_SHRINK:-apeglm}"                   # apeglm | normal | none
INDEPENDENT_FILTERING="${INDEPENDENT_FILTERING:-yes}" # yes | no (default yes, DESeq2's own default; no = pre-1.9 behaviour)

# Named in the report's Methods ("All computational analyses were performed
# on ..."). Recorded at analysis time by the pipeline; override via --export.
# (two-step default: an apostrophe inside "${VAR:-...}" breaks bash's parser)
PLATFORM_NAME="${PLATFORM_NAME:-}"
[ -n "$PLATFORM_NAME" ] || PLATFORM_NAME="VCU's High Performance Research Computing (HPRC) cluster"
export BISR_PLATFORM_NAME="$PLATFORM_NAME"

# How to run: auto | container | conda.
#   container  the prebuilt Apptainer image (bisrde_<version>.sif in BISR_SIF_DIR,
#              or the file named by BISR_SIF). No conda env needed and no MSigDB
#              preflight: the gene sets are baked in. Compute nodes cannot pull,
#              so the image must already exist (pull_container.sh on a login node).
#   conda      the conda environment below, as before.
#   auto       container when apptainer/singularity is on PATH and the image is
#              already present (or BISR_SIF is set), otherwise conda.
RUNTIME="${RUNTIME:-auto}"

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

# --------------------------------------------------------------------------
# Runtime
# --------------------------------------------------------------------------
# Resolve the image the launcher would use, so auto can check it exists.
VERSION="$(sed -n 's/^Version:[[:space:]]*//p' bisrDE/DESCRIPTION 2>/dev/null || true)"
VERSION="${VERSION:-latest}"
SIF_DIR="${BISR_SIF_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/bisrde}"
SIF="${BISR_SIF:-$SIF_DIR/bisrde_$VERSION.sif}"

# Compute nodes have no network, so nothing in this job may pull an image;
# run_analysis.sh says how to fetch it on a login node instead.
export BISR_NO_PULL=1

# An explicit BISR_SIF selects the container whether or not the file exists,
# so a mistyped path fails in the launcher with a message naming it instead
# of silently running conda.
if [ "$RUNTIME" = "auto" ]; then
    if { command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1; } \
       && { [ -n "${BISR_SIF:-}" ] || [ -f "$SIF" ]; }; then
        RUNTIME=container
    else
        RUNTIME=conda
    fi
fi

case "$RUNTIME" in
    container|conda) ;;
    *) echo "ERROR: RUNTIME must be auto, container or conda (got: $RUNTIME)" >&2; exit 2 ;;
esac

if [ "$RUNTIME" = "container" ]; then
    # The launcher does the rest; the image must be there already.
    export BISR_RUNTIME=container
    [ -z "${BISR_SIF_DIR:-}" ] || export BISR_SIF_DIR
    [ -z "${BISR_SIF:-}" ] || export BISR_SIF
    ENV_LINE="$SIF"
    R_LINE="inside the image"
    MSIGDB_STATUS="baked into the image"
else

# Pin the launcher to this decision: its own auto rule would otherwise pick
# the container when a legacy dge_analysis.sif sits next to it or another
# Rscript shadows the env's, after conda has already been activated here.
export BISR_RUNTIME=conda

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

# --------------------------------------------------------------------------
# Preflight: MSigDB cache
# --------------------------------------------------------------------------
# msigdbr downloads its gene-set archive from Zenodo on first use and reuses
# the cached copy afterwards. Compute nodes usually have no outbound internet,
# and a missing cache does not stop the pipeline: Hallmark enrichment times
# out later in the run and the report ends up without it. Check before
# spending the allocation. Resolve the directory exactly as R will
# (tools::R_user_dir honours R_USER_CACHE_DIR); --no-init-file skips the
# project .Rprofile, whose renv/conda banner would otherwise leak into the
# captured path.
MSIGDB_CACHE="$(Rscript --no-init-file -e 'cat(tools::R_user_dir("msigdbr", "cache"))' 2>/dev/null || true)"
if [ -n "$MSIGDB_CACHE" ] && [ -n "$(ls -A "$MSIGDB_CACHE" 2>/dev/null)" ]; then
    MSIGDB_STATUS="$MSIGDB_CACHE"
else
    MSIGDB_STATUS="${MSIGDB_CACHE:-<unresolved>}  (MISSING: Hallmark will fail offline, see the .err log)"
    {
        echo "WARNING: MSigDB cache is missing or empty: ${MSIGDB_CACHE:-<unresolved>}"
        echo "  Compute nodes usually cannot reach zenodo.org, so Hallmark enrichment"
        echo "  will fail with 'Timeout was reached [zenodo.org]' and the report will"
        echo "  have no Hallmark results. Everything else still runs."
        echo "  Warm the cache ONCE on a login node, with the same cache location this"
        echo "  job resolved:"
        [ -z "${R_USER_CACHE_DIR:-}" ] || echo "    export R_USER_CACHE_DIR=$R_USER_CACHE_DIR"
        echo "    cd $PIPELINE_DIR && Rscript warm_msigdb_cache.R"
        echo "  then resubmit."
        echo
    } >&2
fi

ENV_LINE="${CONDA_PREFIX:-<none>}"
R_LINE="$(command -v Rscript) ($(Rscript --no-init-file -e 'cat(R.version.string)' 2>/dev/null))"

fi  # RUNTIME

echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-<interactive>}"
echo "runtime:   $RUNTIME"
echo "env:       $ENV_LINE"
echo "R:         $R_LINE"
echo "outdir:    $OUTDIR"
echo "methods:   gsea-rank=$GSEA_RANK  lfc-shrink=$LFC_SHRINK  independent-filtering=$INDEPENDENT_FILTERING"
echo "platform:  $BISR_PLATFORM_NAME"
echo "msigdb:    $MSIGDB_STATUS"
echo

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
args=(--counts "$COUNTS" --samplesheet "$SAMPLESHEET"
      --outdir "$OUTDIR" --runid "$RUNID"
      --annotation "$ANNOTATION" --id-type "$IDTYPE"
      --fold-change "$FOLD_CHANGE" --padj "$PADJ"
      --gsea-rank "$GSEA_RANK" --lfc-shrink "$LFC_SHRINK"
      --independent-filtering "$INDEPENDENT_FILTERING")
[ -n "$BRS_TICKET" ] && args+=(--brs-ticket "$BRS_TICKET")
[ -n "${ANALYST//[[:space:]]/}" ] && args+=(--analyst "$ANALYST")

bash run_analysis.sh "${args[@]}"

echo
echo "done. Report + logs under: $OUTDIR"
