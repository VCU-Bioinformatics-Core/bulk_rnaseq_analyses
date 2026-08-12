#!/bin/bash
#
# build_container.sh
# Script to build the Singularity/Apptainer container for differential expression analysis
#
# Usage:
#   bash build_container.sh
#
# Requirements:
#   - Apptainer/Singularity installed
#   - Sufficient disk space (~3-4 GB)
#   - May require sudo or --fakeroot flag depending on system setup

set -e  # Exit on error
set -u  # Exit on undefined variable

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "======================================================================"
echo "Differential Expression Analysis - Container Build"
echo "======================================================================"
echo ""

# Check if definition file exists
if [ ! -f "dge_analysis.def" ]; then
    echo -e "${RED}ERROR: dge_analysis.def not found!${NC}"
    echo "Please run this script from the differential_expression directory."
    exit 1
fi

# Check if lockfile exists (Phase 8.1 collapsed dge_renv.lock into renv.lock)
if [ ! -f "renv.lock" ]; then
    echo -e "${RED}ERROR: renv.lock not found!${NC}"
    echo "This file is required for building the container."
    exit 1
fi

# Check for apptainer or singularity
CONTAINER_CMD=""
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
    echo -e "${GREEN}✓ Found: apptainer${NC}"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
    echo -e "${GREEN}✓ Found: singularity${NC}"
else
    echo -e "${RED}ERROR: Neither apptainer nor singularity found!${NC}"
    echo "Please install Apptainer or Singularity to continue."
    exit 1
fi

echo ""
echo "Building container image..."
echo "This will take approximately 45-90 minutes on first build."
echo "The container will be saved as: dge_analysis.sif"
echo ""

# Ask for build method
echo "Select build method:"
echo "  1) Standard build (requires sudo)"
echo "  2) Fakeroot build (no sudo required, if fakeroot is configured)"
echo ""
read -r -p "Enter choice [1/2]: " BUILD_CHOICE

case $BUILD_CHOICE in
    1)
        echo -e "${YELLOW}Building with sudo...${NC}"
        sudo $CONTAINER_CMD build dge_analysis.sif dge_analysis.def
        ;;
    2)
        echo -e "${YELLOW}Building with fakeroot...${NC}"
        $CONTAINER_CMD build --fakeroot dge_analysis.sif dge_analysis.def
        ;;
    *)
        echo -e "${RED}Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

# Check if build was successful
if [ -f "dge_analysis.sif" ]; then
    echo ""
    echo "======================================================================"
    echo -e "${GREEN}✓ Container build successful!${NC}"
    echo "======================================================================"
    echo ""
    echo "Container file: dge_analysis.sif"
    echo "Size: $(du -h dge_analysis.sif | cut -f1)"
    echo ""
    echo "To test the container:"
    echo "  $CONTAINER_CMD exec dge_analysis.sif R --version"
    echo "  $CONTAINER_CMD exec dge_analysis.sif Rscript -e 'library(DESeq2); packageVersion(\"DESeq2\")'"
    echo ""
    echo "To run analysis:"
    echo "  bash run_in_container.sh --counts <file> --samplesheet <file> ..."
    echo ""
else
    echo -e "${RED}ERROR: Container build failed!${NC}"
    exit 1
fi
