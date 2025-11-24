#!/bin/bash
#
# run_in_container.sh  
# Wrapper script to run differential expression analysis inside Singularity/Apptainer container
#
# Usage:
#   bash run_in_container.sh --counts <file> --samplesheet <file> --outdir <dir> --runid <id> --annotation <mouse|human>
#
# Example:
#   bash run_in_container.sh \
#     --counts example_counts.tsv \
#     --samplesheet example_samplesheet.csv \
#     --outdir ./results \
#     --runid test_run \
#     --annotation mouse

set -e
set -u

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "======================================================================"
echo "Differential Expression Analysis - Container Execution"
echo "======================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if container exists
CONTAINER_IMAGE="$SCRIPT_DIR/dge_analysis.sif"
if [ ! -f "$CONTAINER_IMAGE" ]; then
    echo -e "${RED}ERROR: Container image not found: $CONTAINER_IMAGE${NC}"
    echo ""
    echo "Please build the container first:"
    echo "  bash build_container.sh"
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

# Pass all arguments to the R script inside the container
# Bind mount the current directory and common data directories
echo "Running analysis in container..."
echo "Arguments: $@"
echo ""

# Get current directory for binding
CWD=$(pwd)

# Run the analysis
$CONTAINER_CMD exec \
    --bind "$CWD:/workspace" \
    --bind "$SCRIPT_DIR:/scripts" \
    --pwd /workspace \
    "$CONTAINER_IMAGE" \
    Rscript /scripts/de.R "$@"

EXIT_CODE=$?

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
