#!/bin/bash
#
# run_analysis.sh
# Smart wrapper that detects the environment and runs analysis appropriately
#
# Usage:
#   bash run_analysis.sh --counts <file> --samplesheet <file> --outdir <dir> --runid <id> --annotation <mouse|human>
#
# This script will:
#   - On Mac: Use local R with renv
#   - On Linux with container: Use Singularity/Apptainer container
#   - On Linux without container: Use local R with renv

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

# Check if container image exists
CONTAINER_IMAGE="$SCRIPT_DIR/dge_analysis.sif"
HAS_CONTAINER_IMAGE=false
if [ -f "$CONTAINER_IMAGE" ]; then
    HAS_CONTAINER_IMAGE=true
fi

# Decide execution method
EXEC_METHOD=""

if [ "$OS_TYPE" == "Darwin" ]; then
    # macOS - use local R with renv
    echo -e "${BLUE}Platform: macOS - Will use local R with renv${NC}"
    EXEC_METHOD="local"
    
elif [ "$OS_TYPE" == "Linux" ]; then
    # Linux - prefer container if available
    if [ "$HAS_CONTAINER" = true ] && [ "$HAS_CONTAINER_IMAGE" = true ]; then
        echo -e "${BLUE}Platform: Linux - Will use $CONTAINER_CMD container${NC}"
        EXEC_METHOD="container"
    else
        echo -e "${YELLOW}Platform: Linux - Container not available, will use local R with renv${NC}"
        EXEC_METHOD="local"
        
        if [ "$HAS_CONTAINER" = false ]; then
            echo -e "${YELLOW}Note: Install Apptainer/Singularity for containerized execution${NC}"
        elif [ "$HAS_CONTAINER_IMAGE" = false ]; then
            echo -e "${YELLOW}Note: Build container first: bash build_container.sh${NC}"
        fi
    fi
else
    echo -e "${YELLOW}Unknown OS: $OS_TYPE - Will attempt local R with renv${NC}"
    EXEC_METHOD="local"
fi

echo ""
echo "Execution method: $EXEC_METHOD"
echo ""

# Execute based on method
if [ "$EXEC_METHOD" == "container" ]; then
    # Run in container
    echo -e "${GREEN}Launching analysis in container...${NC}"
    exec bash run_in_container.sh "$@"
    
elif [ "$EXEC_METHOD" == "local" ]; then
    # Run locally with renv
    echo -e "${GREEN}Launching analysis with local R...${NC}"
    
    # Check if R is available
    if ! command -v Rscript &> /dev/null; then
        echo -e "${RED}ERROR: Rscript not found!${NC}"
        echo "Please install R to continue."
        exit 1
    fi
    
    # Check if renv is set up
    if [ ! -d "renv" ]; then
        echo -e "${YELLOW}renv not initialized. Initializing now...${NC}"
        echo "This is a one-time setup and may take 30-60 minutes."
        echo ""
        Rscript setup_renv.R
        echo ""
    fi
    
    # Run the analysis
    echo "Running: Rscript de.R $*"
    echo ""
    exec Rscript de.R "$@"
    
else
    echo -e "${RED}ERROR: Could not determine execution method${NC}"
    exit 1
fi
