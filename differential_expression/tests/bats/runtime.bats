#!/usr/bin/env bats
#
# Runtime resolution in the launcher, the container wrapper, the pull script
# and the Slurm template.
#
# Everything runs with BISR_DRY_RUN=1 and stub `uname` / `apptainer` / `Rscript`
# binaries on a throwaway PATH, so no image is pulled and no R is started.
# `uname` is stubbed to print Linux because run_analysis.sh never picks the
# container on macOS.
#
# Run with:  just test-sh     (or: bats tests/bats)

setup() {
    DE="$( cd "$( dirname "$BATS_TEST_FILENAME" )/../.." && pwd -P )"
    export DE
    VERSION="$(sed -n 's/^Version:[[:space:]]*//p' "$DE/bisrDE/DESCRIPTION")"
    [ -n "$VERSION" ]

    TMP="$( cd "$(mktemp -d)" && pwd -P )"
    export TMP

    # Stub PATH: uname always says Linux; nothing else from the host beyond
    # the system directories, so a real apptainer or conda cannot leak in.
    STUBS="$TMP/stubs"; mkdir -p "$STUBS"
    printf '#!/bin/sh\necho Linux\n' > "$STUBS/uname"
    chmod +x "$STUBS/uname"
    export PATH="$STUBS:/usr/bin:/bin"

    # Optional apptainer stub: echoes its argv; `pull <sif> <image>` creates <sif>.
    APPT="$TMP/apptainer-bin"; mkdir -p "$APPT"
    cat > "$APPT/apptainer" <<'EOF'
#!/bin/sh
echo "apptainer $*"
if [ "$1" = pull ]; then : > "$2"; fi
EOF
    chmod +x "$APPT/apptainer"

    export BISR_DRY_RUN=1
    export BISR_SIF_DIR="$TMP/sifs"
    unset BISR_SIF BISR_RUNTIME BISR_NO_PULL BISR_IMAGE BISR_BIND CONDA_PREFIX

    ARGS=(--counts assets/example_counts.tsv --samplesheet assets/example_samplesheet.csv
          --outdir "$TMP/out" --runid rt --annotation mouse)
}

teardown() {
    rm -rf "$TMP"
}

with_apptainer() { export PATH="$APPT:$PATH"; }

launch() {
    ( cd "$DE" && bash run_analysis.sh "${ARGS[@]}" 2>&1 )
}

line() { printf '%s\n' "$output" | grep "^$1:" ; }

@test "BISR_SIF set resolves to container and the command names that image" {
    with_apptainer
    : > "$TMP/custom.sif"
    export BISR_SIF="$TMP/custom.sif"
    run launch
    [ "$status" -eq 0 ]
    [ "$(line runtime)" = "runtime: container" ]
    [[ "$(line command)" == *"$TMP/custom.sif"* ]]
    [[ "$(line command)" == *"run_in_container.sh"* ]]
}

@test "legacy dge_analysis.sif next to the script is still honoured" {
    with_apptainer
    mkdir -p "$TMP/de/bisrDE"
    cp "$DE/run_analysis.sh" "$DE/run_in_container.sh" "$TMP/de/"
    cp "$DE/bisrDE/DESCRIPTION" "$TMP/de/bisrDE/"
    : > "$TMP/de/dge_analysis.sif"
    run bash "$TMP/de/run_analysis.sh" "${ARGS[@]}"
    [ "$status" -eq 0 ]
    [ "$(line runtime)" = "runtime: container" ]
    [[ "$(line command)" == *"$TMP/de/dge_analysis.sif"* ]]
}

@test "active conda env providing Rscript resolves to conda" {
    export CONDA_PREFIX="$TMP/conda"
    mkdir -p "$CONDA_PREFIX/bin"
    printf '#!/bin/sh\nexit 0\n' > "$CONDA_PREFIX/bin/Rscript"
    chmod +x "$CONDA_PREFIX/bin/Rscript"
    export PATH="$CONDA_PREFIX/bin:$PATH"
    run launch
    [ "$status" -eq 0 ]
    [ "$(line runtime)" = "runtime: conda" ]
    [[ "$(line command)" == *"de.R"* ]]
}

@test "missing SIF with BISR_NO_PULL fails and names pull_container.sh" {
    with_apptainer
    export BISR_NO_PULL=1
    run launch
    [ "$status" -ne 0 ]
    [[ "$output" == *"pull_container.sh"* ]]
    [[ "$output" == *"bisrde_$VERSION.sif"* ]]
    [ ! -e "$BISR_SIF_DIR/bisrde_$VERSION.sif" ]
}

@test "missing BISR_SIF with BISR_NO_PULL names BISR_SIF, not the versioned image" {
    with_apptainer
    export BISR_NO_PULL=1 BISR_SIF="$TMP/typo.sif"
    run launch
    [ "$status" -ne 0 ]
    [[ "$output" == *"BISR_SIF points to a file that does not exist: $TMP/typo.sif"* ]]
    [[ "$output" != *"Container image not found"* ]]
}

@test "apptainer present and no SIF resolves to the versioned image without pulling" {
    with_apptainer
    run launch
    [ "$status" -eq 0 ]
    [ "$(line runtime)" = "runtime: container" ]
    [[ "$(line command)" == *"$BISR_SIF_DIR/bisrde_$VERSION.sif"* ]]
    [ ! -e "$BISR_SIF_DIR/bisrde_$VERSION.sif" ]
    [ ! -d "$BISR_SIF_DIR" ]
}

@test "BISR_RUNTIME=renv forces local R even with apptainer present" {
    with_apptainer
    export BISR_RUNTIME=renv
    run launch
    [ "$status" -eq 0 ]
    [ "$(line runtime)" = "runtime: renv" ]
    [[ "$(line command)" != *"run_in_container.sh"* ]]
}

@test "run_in_container.sh binds the input dirs, the outdir and passes --pwd" {
    with_apptainer
    : > "$TMP/image.sif"
    export BISR_SIF="$TMP/image.sif"
    run bash -c "cd '$DE' && bash run_in_container.sh --counts assets/example_counts.tsv \
        --samplesheet=assets/example_samplesheet.csv --outdir '$TMP/out' --runid rt"
    [ "$status" -eq 0 ]
    cmd="$(line command)"
    [[ "$cmd" == *"apptainer exec --bind "* ]]
    [[ "$cmd" == *"$DE/assets"* ]]
    [[ "$cmd" == *"$TMP/out"* ]]
    [[ "$cmd" == *"--pwd $DE "* ]]
    [[ "$cmd" == *"$TMP/image.sif Rscript /opt/bisrde/de.R --counts assets/example_counts.tsv"* ]]
    # A dry run creates nothing.
    [ ! -d "$TMP/out" ]
}

@test "run_in_container.sh binds a symlinked input dir as typed and resolved" {
    with_apptainer
    : > "$TMP/image.sif"
    export BISR_SIF="$TMP/image.sif"
    mkdir -p "$TMP/real data"
    : > "$TMP/real data/c.tsv"
    ln -s "$TMP/real data" "$TMP/link"
    run bash "$DE/run_in_container.sh" --counts "$TMP/link/c.tsv" --samplesheet "$TMP/link/s.csv" --outdir "$TMP/out"
    [ "$status" -eq 0 ]
    binds="$(printf '%s\n' "$output" | sed -n 's/^Binding:  *//p')"
    [[ ",$binds," == *",$TMP/link,"* ]]
    [[ ",$binds," == *",$TMP/real data,"* ]]
    # The argv keeps the path the user typed.
    [[ "$(line command)" == *"--counts $TMP/link/c.tsv"* ]]
}

@test "run_in_container.sh trims whitespace around BISR_BIND entries" {
    with_apptainer
    : > "$TMP/image.sif"
    export BISR_SIF="$TMP/image.sif"
    mkdir -p "$TMP/extra one" "$TMP/extra two"
    export BISR_BIND="$TMP/extra one, $TMP/extra two ,"
    run bash "$DE/run_in_container.sh" --counts "$TMP/c.tsv" --outdir "$TMP/out"
    [ "$status" -eq 0 ]
    binds="$(printf '%s\n' "$output" | sed -n 's/^Binding:  *//p')"
    [[ ",$binds," == *",$TMP/extra one,"* ]]
    [[ ",$binds," == *",$TMP/extra two,"* ]]
    [[ "$binds" != *", "* ]]
}

@test "run_in_container.sh runs the checkout's de.R inside a legacy dge_analysis.sif" {
    with_apptainer
    : > "$TMP/dge_analysis.sif"
    export BISR_SIF="$TMP/dge_analysis.sif"
    run bash "$DE/run_in_container.sh" --counts "$TMP/c.tsv" --outdir "$TMP/out"
    [ "$status" -eq 0 ]
    [[ "$(line command)" == *"$TMP/dge_analysis.sif Rscript $DE/de.R --counts"* ]]
    [[ "$(line command)" != *"/opt/bisrde/de.R"* ]]
    [[ "$output" == *"pull_container.sh"* ]]
}

@test "run_in_container.sh refuses to run without BISR_SIF" {
    with_apptainer
    run bash "$DE/run_in_container.sh" "${ARGS[@]}"
    [ "$status" -ne 0 ]
    [[ "$output" == *"BISR_SIF"* ]]
}

@test "pull_container.sh creates the versioned SIF and is idempotent" {
    with_apptainer
    unset BISR_DRY_RUN
    run bash "$DE/pull_container.sh"
    [ "$status" -eq 0 ]
    [ -f "$BISR_SIF_DIR/bisrde_$VERSION.sif" ]
    [[ "$output" == *"docker://ghcr.io/vcu-bioinformatics-core/bisrde:$VERSION"* ]]
    [[ "$output" == *"export BISR_SIF="* ]]
    run bash "$DE/pull_container.sh"
    [ "$status" -eq 0 ]
    [[ "$output" == *"already present"* ]]
}

@test "pull_container.sh exits 2 without apptainer or singularity" {
    unset BISR_DRY_RUN
    run bash "$DE/pull_container.sh"
    [ "$status" -eq 2 ]
}

@test "submit_slurm.sh RUNTIME=container exports the container knobs and needs no conda" {
    with_apptainer
    : > "$TMP/image.sif"
    mkdir -p "$TMP/pipe"
    printf '#!/bin/sh\nenv | grep "^BISR_" | sort\necho "args: $*"\n' > "$TMP/pipe/run_analysis.sh"
    export PIPELINE_DIR="$TMP/pipe" RUNTIME=container BISR_SIF="$TMP/image.sif"
    unset BISR_DRY_RUN
    run bash "$DE/submit_slurm.sh"
    [ "$status" -eq 0 ]
    [[ "$output" == *"BISR_RUNTIME=container"* ]]
    [[ "$output" == *"BISR_NO_PULL=1"* ]]
    [[ "$output" == *"BISR_SIF=$TMP/image.sif"* ]]
    [[ "$output" == *"msigdb:    baked into the image"* ]]
    [[ "$output" == *"env:       $TMP/image.sif"* ]]
}

@test "submit_slurm.sh RUNTIME=auto picks the container only when the SIF exists" {
    with_apptainer
    mkdir -p "$TMP/pipe/bisrDE"
    cp "$DE/bisrDE/DESCRIPTION" "$TMP/pipe/bisrDE/"
    printf '#!/bin/sh\nenv | grep "^BISR_" | sort\n' > "$TMP/pipe/run_analysis.sh"
    export PIPELINE_DIR="$TMP/pipe"
    unset BISR_DRY_RUN
    mkdir -p "$BISR_SIF_DIR"
    : > "$BISR_SIF_DIR/bisrde_$VERSION.sif"
    run bash "$DE/submit_slurm.sh"
    [ "$status" -eq 0 ]
    [[ "$output" == *"runtime:   container"* ]]
    [[ "$output" == *"BISR_RUNTIME=container"* ]]
    # Without the image, auto falls back to conda, which is absent here.
    rm "$BISR_SIF_DIR/bisrde_$VERSION.sif"
    run bash "$DE/submit_slurm.sh"
    [ "$status" -ne 0 ]
    [[ "$output" == *"conda"* ]]
}

@test "submit_slurm.sh RUNTIME=auto picks the container when BISR_SIF is set but missing" {
    with_apptainer
    mkdir -p "$TMP/pipe"
    printf '#!/bin/sh\nenv | grep "^BISR_" | sort\n' > "$TMP/pipe/run_analysis.sh"
    export PIPELINE_DIR="$TMP/pipe" BISR_SIF="$TMP/typo.sif"
    unset BISR_DRY_RUN
    run bash "$DE/submit_slurm.sh"
    [ "$status" -eq 0 ]
    [[ "$output" == *"runtime:   container"* ]]
    [[ "$output" == *"BISR_SIF=$TMP/typo.sif"* ]]
    [[ "$output" != *"conda activate"* ]]
}

@test "submit_slurm.sh conda branch pins BISR_RUNTIME=conda and BISR_NO_PULL=1" {
    with_apptainer
    mkdir -p "$TMP/pipe" "$TMP/conda/bin" "$TMP/condabin"
    printf '#!/bin/sh\nenv | grep "^BISR_" | sort\n' > "$TMP/pipe/run_analysis.sh"
    printf '#!/bin/sh\nexit 0\n' > "$TMP/conda/bin/Rscript"
    chmod +x "$TMP/conda/bin/Rscript"
    # `conda shell.bash hook` output is eval'd by the template; the stub's hook
    # defines a conda() that activates the fake env.
    printf '#!/bin/sh\necho "conda() { export CONDA_PREFIX=%s; export PATH=%s/bin:\\$PATH; }"\n' \
        "$TMP/conda" "$TMP/conda" > "$TMP/condabin/conda"
    chmod +x "$TMP/condabin/conda"
    export PATH="$TMP/condabin:$PATH" PIPELINE_DIR="$TMP/pipe"
    unset BISR_DRY_RUN
    run bash "$DE/submit_slurm.sh"
    [ "$status" -eq 0 ]
    [[ "$output" == *"runtime:   conda"* ]]
    [[ "$output" == *"BISR_RUNTIME=conda"* ]]
    [[ "$output" == *"BISR_NO_PULL=1"* ]]
}
