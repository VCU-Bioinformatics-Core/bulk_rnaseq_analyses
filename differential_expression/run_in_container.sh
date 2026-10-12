#!/bin/bash
#
# run_in_container.sh
# Wrapper script to run differential expression analysis inside the bisrDE
# Apptainer/Singularity image. run_analysis.sh resolves the image and calls
# this with BISR_SIF set; it can also be called directly.
#
# Usage:
#   BISR_SIF=/path/to/bisrde_<version>.sif \
#   bash run_in_container.sh --counts <file> --samplesheet <file> --outdir <dir> --runid <id> --annotation <mouse|human>
#
# Binds into the container: the directories holding --counts and
# --samplesheet, the --outdir (created first), the current directory, and
# every path listed in BISR_BIND (comma-separated). Each is bound both as
# typed and with symlinks resolved. The driver and package baked into the
# image are what runs (/opt/bisrde/de.R), not the checkout, so the code
# always matches the image version; the one exception is the pre-1.10
# dge_analysis.sif, which has no baked driver and runs the checkout's de.R.
# The host environment passes through unchanged (BISR_PLATFORM_NAME and
# friends reach the run).
#
# BISR_DRY_RUN=1 prints the exact command instead of running it.

set -e
set -u

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

echo "======================================================================"
echo "Differential Expression Analysis - Container Execution"
echo "======================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if container exists
CONTAINER_IMAGE="${BISR_SIF:-}"
if [ -z "$CONTAINER_IMAGE" ]; then
    echo -e "${RED}ERROR: BISR_SIF is not set${NC}"
    echo ""
    echo "Run through the launcher, which resolves the image:"
    echo "  bash run_analysis.sh ..."
    echo "or point at an image explicitly:"
    echo "  BISR_SIF=/path/to/bisrde_<version>.sif bash run_in_container.sh ..."
    echo ""
    exit 1
fi
if [ ! -f "$CONTAINER_IMAGE" ]; then
    echo -e "${RED}ERROR: Container image not found: $CONTAINER_IMAGE${NC}"
    echo ""
    echo "Pull it on a login node first:"
    echo "  bash $SCRIPT_DIR/pull_container.sh"
    echo ""
    exit 1
fi

# Check for apptainer or singularity
CONTAINER_CMD=""
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo -e "${RED}ERROR: Neither apptainer nor singularity found!${NC}"
    exit 1
fi

echo -e "${GREEN}✓ Using: $CONTAINER_CMD${NC}"
echo -e "${GREEN}✓ Container: $CONTAINER_IMAGE${NC}"
echo ""

# Absolute path as the user spelled it, symlinks kept (logical pwd). A path
# that does not exist yet is joined onto $PWD.
abspath() {
    local p="$1"
    case "$p" in
        /*) ;;
        *)  p="$PWD/${p#./}" ;;
    esac
    if [ -d "$p" ]; then
        (cd "$p" && pwd)
    elif [ -d "$(dirname "$p")" ]; then
        printf '%s/%s\n' "$(cd "$(dirname "$p")" && pwd)" "$(basename "$p")"
    else
        printf '%s\n' "$p"
    fi
}

# Symlink-free form of an absolute path: readlink -f where it exists (GNU,
# macOS 12.3+), otherwise cd && pwd -P. Empty when neither can resolve it.
physpath() {
    local p="$1" r
    if r="$(readlink -f -- "$p" 2>/dev/null)" && [ -n "$r" ]; then
        printf '%s\n' "$r"
    elif [ -d "$p" ]; then
        (cd "$p" && pwd -P)
    fi
}

# Read --counts/--samplesheet/--outdir out of the arguments without consuming
# them; the full argv is forwarded to de.R untouched.
COUNTS=""
SAMPLESHEET=""
OUTDIR=""
ARGS=("$@")
NARGS=$#
i=0
while [ "$i" -lt "$NARGS" ]; do
    a="${ARGS[$i]}"
    case "$a" in
        --counts)        i=$((i + 1)); COUNTS="${ARGS[$i]:-}" ;;
        --counts=*)      COUNTS="${a#*=}" ;;
        --samplesheet)   i=$((i + 1)); SAMPLESHEET="${ARGS[$i]:-}" ;;
        --samplesheet=*) SAMPLESHEET="${a#*=}" ;;
        --outdir)        i=$((i + 1)); OUTDIR="${ARGS[$i]:-}" ;;
        --outdir=*)      OUTDIR="${a#*=}" ;;
    esac
    i=$((i + 1))
done

# The output directory must exist before the bind, or Apptainer refuses it.
if [ -n "$OUTDIR" ] && [ -z "${BISR_DRY_RUN:-}" ]; then
    mkdir -p "$OUTDIR"
fi

# Bind list, comma-separated and de-duplicated.
BIND_LIST=""
add_bind() {
    local p="$1"
    [ -n "$p" ] || return 0
    case ",$BIND_LIST," in
        *",$p,"*) return 0 ;;
    esac
    BIND_LIST="${BIND_LIST:+$BIND_LIST,}$p"
}


# Bind a path under the name the user typed and under its symlink-free form.
# Apptainer mounts a bind at the literal source path and does not recreate
# host symlinks, so when /scratch is a link to /gpfs/scratch the argv handed
# to de.R (/scratch/...) only resolves inside the container if that spelling
# is bound too. add_bind drops the duplicate when the two forms agree.
bind_both() {
    local p phys
    p="$(abspath "$1")"
    add_bind "$p"
    phys="$(physpath "$p")"
    [ -z "$phys" ] || add_bind "$phys"
}

[ -z "$COUNTS" ]      || bind_both "$(dirname "$COUNTS")"
[ -z "$SAMPLESHEET" ] || bind_both "$(dirname "$SAMPLESHEET")"
[ -z "$OUTDIR" ]      || bind_both "$OUTDIR"
bind_both "$PWD"
if [ -n "${BISR_BIND:-}" ]; then
    OLD_IFS="$IFS"; IFS=','
    for extra in $BISR_BIND; do
        # Trim whitespace around each entry so "a, b" works.
        extra="${extra#"${extra%%[![:space:]]*}"}"
        extra="${extra%"${extra##*[![:space:]]}"}"
        [ -z "$extra" ] || bind_both "$extra"
    done
    IFS="$OLD_IFS"
fi

# The driver baked into the image. The pre-1.10 dge_analysis.sif installed
# bisrDE but not de.R, and its wrapper ran the checkout's driver through a
# bind; keep doing that for it.
DRIVER=/opt/bisrde/de.R
if [ "$(basename "$CONTAINER_IMAGE")" = dge_analysis.sif ]; then
    DRIVER="$SCRIPT_DIR/de.R"
    bind_both "$SCRIPT_DIR"
    echo "NOTE: legacy dge_analysis.sif; running the checkout's de.R against the packages in that image."
    echo "      Pull the versioned image with: bash $SCRIPT_DIR/pull_container.sh"
    echo ""
fi

CMD=("$CONTAINER_CMD" exec --bind "$BIND_LIST" --pwd "$PWD" "$CONTAINER_IMAGE"
     Rscript "$DRIVER" "$@")

echo "Binding:   $BIND_LIST"
echo "Arguments: $*"
echo ""

if [ -n "${BISR_DRY_RUN:-}" ]; then
    echo "runtime: container"
    echo "command: ${CMD[*]}"
    exit 0
fi

echo "Running analysis in container..."
echo ""

EXIT_CODE=0
"${CMD[@]}" || EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "======================================================================"
    echo -e "${GREEN}✓ Analysis complete!${NC}"
    echo "======================================================================"
else
    echo "======================================================================"
    echo -e "${RED}✗ Analysis failed with exit code: $EXIT_CODE${NC}"
    echo "======================================================================"
fi

exit $EXIT_CODE
