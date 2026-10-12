#!/bin/bash
#
# run_analysis.sh
# Smart wrapper that detects the environment and runs analysis appropriately
#
# Usage:
#   bash run_analysis.sh --counts <file> --samplesheet <file> --outdir <dir> --runid <id> --annotation <mouse|human>
#
# Runtime selection (BISR_RUNTIME=auto, the default):
#   - macOS: conda env if one is active, otherwise local R with renv.
#     Apptainer does not run on macOS, so the container is never chosen.
#   - Linux, in order:
#       BISR_SIF set                                  -> container
#       legacy dge_analysis.sif next to this script   -> container
#       active conda env providing Rscript            -> conda
#       apptainer or singularity on PATH              -> container
#                                                        (pulls the versioned image on first use)
#       otherwise                                     -> local R with renv
#
# Environment knobs:
#   BISR_RUNTIME   auto | container | conda | renv   (default auto)
#   BISR_SIF       explicit .sif path; wins over everything
#   BISR_SIF_DIR   where pulled images live (default ${XDG_CACHE_HOME:-$HOME/.cache}/bisrde);
#                  the image inside it is bisrde_<version>.sif
#   BISR_IMAGE     image reference to pull (default docker://ghcr.io/vcu-bioinformatics-core/bisrde:<version>)
#   BISR_NO_PULL   never pull (compute nodes); errors with the pull command instead
#   BISR_BIND      extra comma-separated host paths to bind into the container
#   BISR_DRY_RUN   print "runtime: ..." and "command: ..." and exit without running
#
# BISR_NO_PULL and BISR_DRY_RUN are on when set to any non-empty value (1 by
# convention); unset them to turn them off.
# <version> is the Version field of bisrDE/DESCRIPTION.
#
# The launcher changes into this directory before running anything, in every
# runtime, so relative --counts/--samplesheet/--outdir paths resolve against
# the differential_expression directory, not the caller's.

set -e
set -u

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo "======================================================================"
echo "Differential Expression Analysis - Launcher"
echo "======================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Detect operating system
OS_TYPE=$(uname -s)
echo "Detected OS: $OS_TYPE"

# Pipeline version, which names the image tag and the SIF file.
VERSION="$(sed -n 's/^Version:[[:space:]]*//p' "$SCRIPT_DIR/bisrDE/DESCRIPTION" 2>/dev/null || true)"
VERSION="${VERSION:-latest}"

BISR_RUNTIME="${BISR_RUNTIME:-auto}"
BISR_SIF_DIR="${BISR_SIF_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/bisrde}"
BISR_IMAGE="${BISR_IMAGE:-docker://ghcr.io/vcu-bioinformatics-core/bisrde:$VERSION}"
LEGACY_SIF="$SCRIPT_DIR/dge_analysis.sif"
VERSIONED_SIF="$BISR_SIF_DIR/bisrde_$VERSION.sif"

# Check for Singularity/Apptainer
HAS_CONTAINER=false
CONTAINER_CMD=""
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
    HAS_CONTAINER=true
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
    HAS_CONTAINER=true
fi

# An active conda env counts only when it is the one providing Rscript.
CONDA_ACTIVE=false
if [ -n "${CONDA_PREFIX:-}" ] && \
   [ "$(command -v Rscript 2>/dev/null || true)" = "${CONDA_PREFIX}/bin/Rscript" ]; then
    CONDA_ACTIVE=true
fi

# Decide execution method: container | conda | renv
EXEC_METHOD=""

case "$BISR_RUNTIME" in
    container|conda|renv)
        echo -e "${BLUE}Runtime forced by BISR_RUNTIME=$BISR_RUNTIME${NC}"
        EXEC_METHOD="$BISR_RUNTIME"
        ;;
    auto)
        if [ "$OS_TYPE" == "Darwin" ]; then
            # macOS - Apptainer does not run here, so local R only
            if [ "$CONDA_ACTIVE" = true ]; then
                echo -e "${BLUE}Platform: macOS - Will use the active conda env${NC}"
                EXEC_METHOD="conda"
            else
                echo -e "${BLUE}Platform: macOS - Will use local R with renv${NC}"
                EXEC_METHOD="renv"
            fi

        elif [ "$OS_TYPE" == "Linux" ]; then
            if [ -n "${BISR_SIF:-}" ]; then
                echo -e "${BLUE}Platform: Linux - Will use the container in BISR_SIF${NC}"
                EXEC_METHOD="container"
            elif [ -f "$LEGACY_SIF" ]; then
                echo -e "${BLUE}Platform: Linux - Will use the legacy container dge_analysis.sif${NC}"
                EXEC_METHOD="container"
            elif [ "$CONDA_ACTIVE" = true ]; then
                echo -e "${BLUE}Platform: Linux - Will use the active conda env${NC}"
                EXEC_METHOD="conda"
            elif [ "$HAS_CONTAINER" = true ]; then
                echo -e "${BLUE}Platform: Linux - Will use $CONTAINER_CMD container${NC}"
                EXEC_METHOD="container"
            else
                echo -e "${YELLOW}Platform: Linux - No conda env and no Apptainer/Singularity, will use local R with renv${NC}"
                echo -e "${YELLOW}Note: Install Apptainer/Singularity, or activate the bisrde conda env, to skip the renv setup${NC}"
                EXEC_METHOD="renv"
            fi
        else
            echo -e "${YELLOW}Unknown OS: $OS_TYPE - Will attempt local R${NC}"
            if [ "$CONDA_ACTIVE" = true ]; then EXEC_METHOD="conda"; else EXEC_METHOD="renv"; fi
        fi
        ;;
    *)
        echo -e "${RED}ERROR: BISR_RUNTIME must be auto, container, conda or renv (got: $BISR_RUNTIME)${NC}"
        exit 1
        ;;
esac

echo ""
echo "Execution method: $EXEC_METHOD"
echo ""

# Execute based on method
if [ "$EXEC_METHOD" == "container" ]; then

    if [ "$HAS_CONTAINER" = false ]; then
        echo -e "${RED}ERROR: Neither apptainer nor singularity found on PATH${NC}"
        echo "Install Apptainer, or set BISR_RUNTIME=conda / BISR_RUNTIME=renv to run without it."
        exit 1
    fi

    # Resolve the image: explicit path, then the legacy file, then the versioned one.
    if [ -n "${BISR_SIF:-}" ]; then
        SIF="$BISR_SIF"
    elif [ -f "$LEGACY_SIF" ]; then
        SIF="$LEGACY_SIF"
    else
        SIF="$VERSIONED_SIF"
    fi

    if [ ! -f "$SIF" ]; then
        # An explicit BISR_SIF is never pulled into existence, so say so first:
        # the pull hint below would create a different file.
        if [ -n "${BISR_SIF:-}" ]; then
            echo -e "${RED}ERROR: BISR_SIF points to a file that does not exist: $SIF${NC}"
            echo "Fix the path, or unset BISR_SIF and use the versioned image instead:"
            echo ""
            echo "  BISR_SIF_DIR=$BISR_SIF_DIR bash $SCRIPT_DIR/pull_container.sh $VERSION"
            echo ""
            exit 1
        fi
        if [ -n "${BISR_NO_PULL:-}" ]; then
            echo -e "${RED}ERROR: Container image not found: $SIF${NC}"
            echo "BISR_NO_PULL is set, so it will not be downloaded here (compute nodes have no network)."
            echo "Pull it once on a login node, then resubmit:"
            echo ""
            echo "  BISR_SIF_DIR=$BISR_SIF_DIR bash $SCRIPT_DIR/pull_container.sh $VERSION"
            echo ""
            exit 1
        fi
    fi

    if [ -n "${BISR_DRY_RUN:-}" ]; then
        echo "runtime: container"
        echo "command: BISR_SIF=$SIF bash $SCRIPT_DIR/run_in_container.sh $*"
        exit 0
    fi

    if [ ! -f "$SIF" ]; then
        echo -e "${YELLOW}Container image not present yet: $SIF${NC}"
        echo "Pulling $BISR_IMAGE"
        echo "The first pull downloads about 3 to 4 GB and converts it to a SIF file;"
        echo "this takes several minutes and only happens once per version."
        echo ""
        mkdir -p "$BISR_SIF_DIR"
        # Keep the OCI layer cache next to the image, not in $HOME, which is small on most clusters.
        export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$BISR_SIF_DIR/.cache}"
        export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$APPTAINER_CACHEDIR}"
        mkdir -p "$APPTAINER_CACHEDIR"
        if ! "$CONTAINER_CMD" pull "$SIF" "$BISR_IMAGE"; then
            echo -e "${RED}ERROR: pull failed for $BISR_IMAGE${NC}"
            echo "Check network access and retry, or pull explicitly with: bash $SCRIPT_DIR/pull_container.sh $VERSION"
            exit 1
        fi
        echo ""
    fi

    echo -e "${GREEN}Launching analysis in container...${NC}"
    export BISR_SIF="$SIF"
    exec bash run_in_container.sh "$@"

elif [ "$EXEC_METHOD" == "conda" ] || [ "$EXEC_METHOD" == "renv" ]; then

    if [ -n "${BISR_DRY_RUN:-}" ]; then
        echo "runtime: $EXEC_METHOD"
        echo "command: Rscript $SCRIPT_DIR/de.R $*"
        exit 0
    fi

    echo -e "${GREEN}Launching analysis with local R...${NC}"

    # Check if R is available
    if ! command -v Rscript &> /dev/null; then
        echo -e "${RED}ERROR: Rscript not found!${NC}"
        echo "Please install R to continue."
        exit 1
    fi

    # If R is coming from an activated conda environment, that environment owns
    # the package library. The repo ships a .Rprofile + renv/activate.R, so
    # without this renv autoloads and shadows the conda library — the classic
    # "works for me" breakage on HPC. Disabling the autoloader makes
    # `conda activate bisrde && bash run_analysis.sh ...` work as expected.
    if [ "$EXEC_METHOD" == "conda" ]; then
        if [ -z "${CONDA_PREFIX:-}" ]; then
            echo -e "${RED}ERROR: BISR_RUNTIME=conda but no conda environment is active${NC}"
            exit 1
        fi
        if [ "$(command -v Rscript)" != "${CONDA_PREFIX}/bin/Rscript" ]; then
            echo -e "${YELLOW}NOTE: \$CONDA_PREFIX is set but Rscript is not from that env.${NC}"
            echo "      Using: $(command -v Rscript)"
        else
            echo -e "${GREEN}Using R from conda env: ${CONDA_PREFIX}${NC}"
        fi
        export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
    else
        if [ -n "${CONDA_PREFIX:-}" ]; then
            echo -e "${YELLOW}NOTE: \$CONDA_PREFIX is set but Rscript is not from that env.${NC}"
            echo "      Using: $(command -v Rscript)"
        fi
        # Check if renv is set up (skipped when a conda env is providing R).
        if [ -z "${CONDA_PREFIX:-}" ] && [ ! -d "renv" ]; then
            echo -e "${YELLOW}renv not initialized. Initializing now...${NC}"
            echo "This is a one-time setup and may take 30-60 minutes."
            echo ""
            Rscript setup_renv.R
            echo ""
        fi
    fi

    # Run the analysis
    echo "Running: Rscript de.R $*"
    echo ""
    exec Rscript de.R "$@"

else
    echo -e "${RED}ERROR: Could not determine execution method${NC}"
    exit 1
fi
