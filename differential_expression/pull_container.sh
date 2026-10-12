#!/usr/bin/env bash
#
# pull_container.sh [version]
#
# Pull the prebuilt bisrDE image from GHCR into BISR_SIF_DIR as
# bisrde_<version>.sif, ready for run_analysis.sh and submit_slurm.sh.
# Run it once per version on a login node (compute nodes have no network).
#
#   version        image tag to pull; default is the Version field of bisrDE/DESCRIPTION
#   BISR_SIF_DIR   destination directory (default ${XDG_CACHE_HOME:-$HOME/.cache}/bisrde);
#                  point it at a shared, group-readable directory so the team pulls once
#   BISR_IMAGE     image reference (default docker://ghcr.io/vcu-bioinformatics-core/bisrde:<version>)
#
# APPTAINER_CACHEDIR is set to <BISR_SIF_DIR>/.cache unless already set, so the
# OCI layer cache lands next to the image rather than filling $HOME.
#
# Exit codes: 0 pulled or already present, 1 pull failed, 2 no apptainer/singularity.

set -euo pipefail

usage() {
    echo "usage: bash pull_container.sh [version]    (default version: bisrDE/DESCRIPTION)"
}

case "${1:-}" in
    -h|--help) usage; exit 0 ;;
esac

if [ "$(uname -s)" = "Darwin" ]; then
    echo "pull_container.sh: Apptainer does not run on macOS; pull on a Linux login node instead." >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VERSION="${1:-}"
if [ -z "$VERSION" ]; then
    VERSION="$(sed -n 's/^Version:[[:space:]]*//p' "$SCRIPT_DIR/bisrDE/DESCRIPTION" 2>/dev/null || true)"
    VERSION="${VERSION:-latest}"
fi

BISR_SIF_DIR="${BISR_SIF_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/bisrde}"
BISR_IMAGE="${BISR_IMAGE:-docker://ghcr.io/vcu-bioinformatics-core/bisrde:$VERSION}"
SIF="$BISR_SIF_DIR/bisrde_$VERSION.sif"

job_line() {
    echo ""
    echo "In a job script or shell, point the launcher at it with:"
    echo "  export BISR_SIF=$SIF"
    echo "or, when the whole team shares this directory:"
    echo "  export BISR_SIF_DIR=$BISR_SIF_DIR"
}

if [ -f "$SIF" ]; then
    echo "already present: $SIF"
    job_line
    exit 0
fi

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_CMD=apptainer
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_CMD=singularity
else
    echo "ERROR: neither apptainer nor singularity is on PATH (try: module load apptainer)" >&2
    exit 2
fi

mkdir -p "$BISR_SIF_DIR"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$BISR_SIF_DIR/.cache}"
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$APPTAINER_CACHEDIR}"
mkdir -p "$APPTAINER_CACHEDIR"

echo "pulling:   $BISR_IMAGE"
echo "into:      $SIF"
echo "cache:     $APPTAINER_CACHEDIR"
echo "The download is about 3 to 4 GB and the SIF conversion takes several minutes."
echo ""

if ! "$CONTAINER_CMD" pull "$SIF" "$BISR_IMAGE"; then
    echo "" >&2
    echo "ERROR: pull failed for $BISR_IMAGE" >&2
    rm -f "$SIF"
    exit 1
fi

echo ""
echo "pulled:    $SIF"
job_line
exit 0
