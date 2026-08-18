# bisrDE — one discoverable surface over this repo's entry points.
#
#   just            list every recipe
#   just demo       end-to-end smoke run on the bundled example data
#   just test       R + Go + shell test suites
#
# External tools are optional: every recipe that needs one checks first and
# tells you how to install it instead of failing obscurely. Nothing here is
# required to run the pipeline — `run_analysis.sh` remains the contract.

set shell := ["bash", "-eu", "-o", "pipefail", "-c"]

de      := justfile_directory() / "differential_expression"
rscript := env_var_or_default("BISR_RSCRIPT", "Rscript")
outdir  := env_var_or_default("BISR_OUTDIR", justfile_directory() / "results")

# List all recipes (default).
default:
    @just --list --unsorted

# Fail with an install hint when a tool is missing.
_need tool hint:
    @command -v {{tool}} >/dev/null 2>&1 || { \
        printf '\033[33m!\033[0m %s is not installed.\n  install: %s\n' '{{tool}}' '{{hint}}'; exit 127; }

# ---------------------------------------------------------------- run --------

# Guided launcher (bash session or Go TUI chooser).
run:
    cd {{de}} && bash run_interactive.sh

# Install the `bisrde` wrapper so recipes run from ANY directory.
# (`just` alone only works inside the repo — it searches upward for a justfile.)
install-cli:
    #!/usr/bin/env bash
    set -euo pipefail
    mkdir -p "$HOME/.local/bin"
    ln -sf "{{justfile_directory()}}/bin/bisrde" "$HOME/.local/bin/bisrde"
    printf '\033[32m✓\033[0m installed: %s\n' "$HOME/.local/bin/bisrde"
    case ":$PATH:" in
      *":$HOME/.local/bin:"*) printf '  now run `bisrde tui` from anywhere.\n' ;;
      *) printf '\033[33m!\033[0m add to PATH: export PATH="$HOME/.local/bin:$PATH"\n' ;;
    esac

# Build and run the Go TUI directly.
tui:
    cd {{de}}/tui && go build -o bisrde-tui . && BISR_PROJECT_DIR={{de}} ./bisrde-tui

# Print the assembled run_analysis.sh command without running it.
print-cmd:
    cd {{de}} && bash run_interactive.sh --print-cmd

# End-to-end smoke run on the bundled example data.
demo:
    cd {{de}} && bash run_analysis.sh \
        --counts assets/example_counts.tsv \
        --samplesheet assets/example_samplesheet.csv \
        --outdir {{outdir}} --runid demo_$(date +%Y%m%d_%H%M%S) \
        --annotation mouse --id-type ensembl

# --------------------------------------------------------------- test --------

# Everything: R package, Go TUI, shell launchers.
test: test-r test-go test-sh

# bisrDE testthat suite.
test-r:
    cd {{de}} && {{rscript}} -e 'devtools::test("bisrDE")'

# Go TUI unit tests.
test-go:
    cd {{de}}/tui && go test ./...

# Shell launcher tests (bats) — asserts bash/Go TUI flag parity.
test-sh: (_need "bats" "brew install bats-core")
    cd {{de}} && bats tests/bats

# Full R CMD check of the package.
check:
    cd {{de}} && {{rscript}} -e 'devtools::check("bisrDE")'

# Regenerate roxygen docs + NAMESPACE.
doc:
    cd {{de}} && {{rscript}} -e 'devtools::document("bisrDE")'

# --------------------------------------------------------------- lint --------

# Static analysis for every shell script (catches real quoting/unset bugs).
lint: (_need "shellcheck" "brew install shellcheck")
    cd {{de}} && shellcheck -x *.sh

# R static analysis.
lint-r:
    cd {{de}} && {{rscript}} -e 'lintr::lint_package("bisrDE")'

# Go vet.
lint-go:
    cd {{de}}/tui && go vet ./...

# ---------------------------------------------------------- inspect ----------

# Open a results CSV / counts matrix in visidata (works over SSH on the HPC).
inspect file: (_need "vd" "brew install visidata  # or: pipx install visidata")
    vd {{file}}

# Render a run's session log to colour-preserved HTML next to it.
log-html log: (_need "aha" "brew install aha")
    #!/usr/bin/env bash
    set -euo pipefail
    src='{{log}}'
    out="${src%.log}.html"
    aha --black --title "bisrDE session log" -f "$src" > "$out"
    printf '\033[32m✓\033[0m wrote %s\n' "$out"

# Disk usage for a directory — HPC scratch quotas are a classic silent failure.
disk dir=".": (_need "dust" "brew install dust  # or: brew install ncdu")
    dust {{dir}}

# ------------------------------------------------------------ artifacts ------

# Screenshot a rendered HTML report (desktop + mobile) for docs / visual diffing.
shot report: (_need "pageres" "npm install -g pageres-cli")
    pageres {{report}} 1366x768 375x812 --filename='report-<%= size %>'

# Losslessly shrink the run's PNGs before the self-contained HTML embeds them.
optimize figdir: (_need "optimizt" "npm install -g @funboxteam/optimizt")
    optimizt --lossless {{figdir}}

# Re-record the terminal demo deterministically from the .tape script.
vhs: (_need "vhs" "brew install vhs")
    cd {{de}} && vhs ../docs/demo.tape

# ------------------------------------------------------------- release -------

# Bump the version everywhere (DESCRIPTION, README, launcher, Go TUI banner).
bump version:
    #!/usr/bin/env bash
    set -euo pipefail
    cd {{de}}
    cur=$(awk '/^Version:/{print $2}' bisrDE/DESCRIPTION)
    echo "bumping $cur -> {{version}}"
    sed -i '' "s/^Version: .*/Version: {{version}}/" bisrDE/DESCRIPTION
    sed -i '' "s/v$cur/v{{version}}/g" README.md run_interactive.sh tui/main.go
    grep -rn "{{version}}" bisrDE/DESCRIPTION README.md run_interactive.sh tui/main.go | head
    echo "now update CHANGELOG.md, then: just test && git commit && git tag v{{version}}"

# ---------------------------------------------------------- container --------

# Build the Apptainer image (needs a Linux host with apptainer/singularity).
container:
    cd {{de}} && bash build_container.sh

# Pre-download the MSigDB gene sets so Hallmark enrichment works OFFLINE.
# msigdbr fetches its archive from Zenodo on first use; HPC compute nodes
# usually have no outbound internet, so run this once on a login node.
msigdb-cache:
    #!/usr/bin/env bash
    set -euo pipefail
    cd {{de}}
    # Under conda the env owns the library, so renv's autoloader must be off.
    # Outside conda, renv IS the library — leave it alone. (Same rule as
    # run_analysis.sh.)
    [ -n "${CONDA_PREFIX:-}" ] && export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
    {{rscript}} -e '
      cache <- tools::R_user_dir("msigdbr", "cache")
      cat("cache dir:", cache, "\n")
      for (sp in c("Homo sapiens", "Mus musculus")) {
        n <- nrow(msigdbr::msigdbr(species = sp, collection = "H"))
        cat(sprintf("  %-14s H collection rows: %d\n", sp, n))
      }
      cat("cached files:\n"); print(list.files(cache))
    '
    printf '\033[32m✓\033[0m compute nodes will now read this cache instead of downloading.\n'
    printf '  If $HOME is not shared with the compute nodes, set R_USER_CACHE_DIR\n'
    printf '  to a shared path in your job script and re-run this recipe.\n'

# Create the conda environment from environment.yml (Bioconductor 3.18 / R 4.3).
conda-env:
    #!/usr/bin/env bash
    set -euo pipefail
    cd {{de}}
    if command -v micromamba >/dev/null 2>&1; then solver=micromamba
    elif command -v mamba >/dev/null 2>&1; then solver=mamba
    elif command -v conda >/dev/null 2>&1; then solver=conda
    else printf '! no conda/mamba/micromamba found.\n  install: brew install micromamba\n'; exit 127; fi
    printf 'creating env with %s (channels: conda-forge, bioconda — no defaults)\n' "$solver"
    "$solver" create -n bisrde -f environment.yml --override-channels -c conda-forge -c bioconda -y
    printf '\033[32m✓\033[0m then: conda activate bisrde && R -e '"'"'remotes::install_local("bisrDE")'"'"'\n'

# Diff two runs (e.g. renv reference vs conda candidate) before trusting a new
# environment. Exit 0 = PASS, 1 = WARN, 2 = FAIL, so it can gate a deployment.
compare a b out="":
    #!/usr/bin/env bash
    set -euo pipefail
    cd {{de}}
    extra=""
    [ -n "{{out}}" ] && extra="--out {{out}}"
    # shellcheck disable=SC2086
    {{rscript}} compare_runs.R --a "{{a}}" --b "{{b}}" $extra

# Dry-run the conda solve for a target platform without installing anything.
conda-solve platform="linux-64":
    cd {{de}} && CONDA_OVERRIDE_GLIBC=2.28 CONDA_OVERRIDE_LINUX=5.4 \
        micromamba create -n bisrde_solvecheck --dry-run --platform {{platform}} \
        --override-channels -c conda-forge -c bioconda -f environment.yml -y
