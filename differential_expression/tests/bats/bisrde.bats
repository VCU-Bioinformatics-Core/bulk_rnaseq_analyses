#!/usr/bin/env bats
#
# The `bisrde` command (bin/bisrde) and its modulefile (modulefiles/bisrde).
#
# Everything runs with BISR_DRY_RUN=1 and stub `uname` / `apptainer` / `sbatch`
# binaries on a throwaway PATH, so no image is pulled, no job is submitted and
# no R is started. `uname` is stubbed to print Linux because run_analysis.sh
# never picks the container on macOS. The command is always started from a
# scratch directory outside the clone, as a lab member would.
#
# Run with:  just test-sh     (or: bats tests/bats)

setup() {
    DE="$( cd "$( dirname "$BATS_TEST_FILENAME" )/../.." && pwd -P )"
    export DE
    BISRDE="$DE/bin/bisrde"
    VERSION="$(sed -n 's/^Version:[[:space:]]*//p' "$DE/bisrDE/DESCRIPTION")"
    [ -n "$VERSION" ]

    TMP="$( cd "$(mktemp -d)" && pwd -P )"
    export TMP
    WORK="$TMP/work"; mkdir -p "$WORK"

    STUBS="$TMP/stubs"; mkdir -p "$STUBS"
    printf '#!/bin/sh\necho Linux\n' > "$STUBS/uname"
    # sbatch records one argument per line, so the tests can compare exactly.
    cat > "$STUBS/sbatch" <<'EOF2'
#!/bin/sh
for a in "$@"; do echo "sbatch-arg: $a"; done
EOF2
    chmod +x "$STUBS/uname" "$STUBS/sbatch"
    export PATH="$STUBS:/usr/bin:/bin"

    APPT="$TMP/apptainer-bin"; mkdir -p "$APPT"
    cat > "$APPT/apptainer" <<'EOF2'
#!/bin/sh
echo "apptainer $*"
if [ "$1" = pull ]; then : > "$2"; fi
EOF2
    chmod +x "$APPT/apptainer"

    export BISR_DRY_RUN=1
    export BISR_SIF_DIR="$TMP/sifs"
    unset BISR_SIF BISR_RUNTIME BISR_NO_PULL BISR_IMAGE BISR_BIND CONDA_PREFIX
}

teardown() {
    rm -rf "$TMP"
}

with_apptainer() { export PATH="$APPT:$PATH"; }

# Run bisrde from the scratch directory.
bisrde() { ( cd "$WORK" && "$BISRDE" "$@" 2>&1 ); }

line() { printf '%s\n' "$output" | grep "^$1:" ; }

sbatch_args() { printf '%s\n' "$output" | sed -n 's/^sbatch-arg: //p'; }

# ---------------------------------------------------------------- dispatch

@test "no arguments prints usage to stderr and exits 2" {
    run bash -c "cd '$WORK' && '$BISRDE' 2>&1 >/dev/null"
    [ "$status" -eq 2 ]
    [[ "$output" == *"usage: bisrde"* ]]
    run bash -c "cd '$WORK' && '$BISRDE' 2>/dev/null"
    [ "$status" -eq 2 ]
    [ -z "$output" ]
}

@test "help, -h and --help print usage to stdout and exit 0" {
    for w in help -h --help; do
        run bash -c "cd '$WORK' && '$BISRDE' $w 2>/dev/null"
        [ "$status" -eq 0 ]
        [[ "$output" == *"usage: bisrde"* ]]
    done
}

@test "unknown subcommand exits 2 with an error and the usage" {
    run bisrde frobnicate
    [ "$status" -eq 2 ]
    [[ "$output" == *"unknown command: frobnicate"* ]]
    [[ "$output" == *"usage: bisrde"* ]]
}

@test "a symlink to the command still finds the clone" {
    ln -s "$BISRDE" "$TMP/linked-bisrde"
    run bash -c "cd '$WORK' && '$TMP/linked-bisrde' version"
    [ "$status" -eq 0 ]
    [[ "$output" == "bisrDE $VERSION"* ]]
}

# ----------------------------------------------------------------- version

@test "version prints the DESCRIPTION version and a missing image" {
    run bisrde version
    [ "$status" -eq 0 ]
    [ "$(printf '%s\n' "$output" | sed -n 1p)" = "bisrDE $VERSION" ]
    [ "$(printf '%s\n' "$output" | sed -n 2p)" = "image: $BISR_SIF_DIR/bisrde_$VERSION.sif (missing)" ]
}

@test "version reports the image as present once the file exists" {
    mkdir -p "$BISR_SIF_DIR"
    : > "$BISR_SIF_DIR/bisrde_$VERSION.sif"
    run bisrde version
    [ "$status" -eq 0 ]
    [ "$(printf '%s\n' "$output" | sed -n 2p)" = "image: $BISR_SIF_DIR/bisrde_$VERSION.sif (present)" ]
}

@test "version prefers BISR_SIF over BISR_SIF_DIR" {
    : > "$TMP/custom.sif"
    export BISR_SIF="$TMP/custom.sif"
    run bisrde version
    [ "$status" -eq 0 ]
    [ "$(printf '%s\n' "$output" | sed -n 2p)" = "image: $TMP/custom.sif (present)" ]
}

@test "version without BISR_SIF_DIR falls back to the cache directory" {
    unset BISR_SIF_DIR
    export XDG_CACHE_HOME="$TMP/xdg"
    run bisrde version
    [ "$status" -eq 0 ]
    [ "$(printf '%s\n' "$output" | sed -n 2p)" = "image: $TMP/xdg/bisrde/bisrde_$VERSION.sif (missing)" ]
}

# --------------------------------------------------------------------- run

@test "run resolves relative paths against the caller's directory" {
    with_apptainer
    run bisrde run --counts rel/c.tsv --samplesheet=ss.csv -o out2 --runid r1 --annotation mouse
    [ "$status" -eq 0 ]
    cmd="$(line command)"
    [[ "$cmd" == *"run_in_container.sh --counts $WORK/rel/c.tsv --samplesheet=$WORK/ss.csv -o $WORK/out2 --runid r1 --annotation mouse" ]]
    [[ "$cmd" != *"--outdir"* ]]
}

@test "run without an outdir appends --outdir <cwd>/output after the other flags" {
    with_apptainer
    run bisrde run -c ./c.tsv -s ss.csv --runid r1
    [ "$status" -eq 0 ]
    [[ "$(line command)" == *"run_in_container.sh -c $WORK/c.tsv -s $WORK/ss.csv --runid r1 --outdir $WORK/output" ]]
}

@test "run leaves absolute paths alone and resolves --outdir=" {
    with_apptainer
    run bisrde run --counts /data/c.tsv --samplesheet /data/ss.csv --outdir=res --runid r1
    [ "$status" -eq 0 ]
    [[ "$(line command)" == *"run_in_container.sh --counts /data/c.tsv --samplesheet /data/ss.csv --outdir=$WORK/res --runid r1" ]]
}

@test "run keeps unrelated flags and their order untouched" {
    with_apptainer
    run bisrde run --fold-change 2 --counts c.tsv --padj 0.01 --analyst "A B" --outdir o
    [ "$status" -eq 0 ]
    [[ "$(line command)" == *"--fold-change 2 --counts $WORK/c.tsv --padj 0.01 --analyst A B --outdir $WORK/o" ]]
}

# ------------------------------------------------------------------ submit

@test "submit passes sbatch options, one --export list and the template last" {
    run bisrde submit --time=01:00:00 COUNTS=rel/c.tsv RUNID=r1
    [ "$status" -eq 0 ]
    [ "$(sbatch_args)" = "--time=01:00:00
--export=ALL,PIPELINE_DIR=$DE,COUNTS=$WORK/rel/c.tsv,RUNID=r1
$DE/submit_slurm.sh" ]
}

@test "submit makes SAMPLESHEET and OUTDIR absolute and keeps absolute values" {
    run bisrde submit SAMPLESHEET=./ss.csv OUTDIR=res COUNTS=/data/c.tsv ANNOTATION=human
    [ "$status" -eq 0 ]
    [ "$(sbatch_args | sed -n 1p)" = "--export=ALL,PIPELINE_DIR=$DE,SAMPLESHEET=$WORK/ss.csv,OUTDIR=$WORK/res,COUNTS=/data/c.tsv,ANNOTATION=human" ]
}

@test "submit without settings still names the clone" {
    run bisrde submit
    [ "$status" -eq 0 ]
    [ "$(sbatch_args)" = "--export=ALL,PIPELINE_DIR=$DE
$DE/submit_slurm.sh" ]
}

@test "submit keeps a value with a space as one argument" {
    run bisrde submit "COUNTS=my dir/c.tsv"
    [ "$status" -eq 0 ]
    [ "$(sbatch_args | sed -n 1p)" = "--export=ALL,PIPELINE_DIR=$DE,COUNTS=$WORK/my dir/c.tsv" ]
}

@test "run leaves an empty value empty" {
    with_apptainer
    run bisrde run --counts "" --outdir=
    [ "$status" -eq 0 ]
    [[ "$(line command)" == *"run_in_container.sh --counts  --outdir=" ]]
}

@test "submit rejects a value with a comma" {
    run bisrde submit RUNID=a,b
    [ "$status" -eq 2 ]
    [[ "$output" == *"comma"* ]]
    [[ "$output" != *"sbatch-arg"* ]]
}

@test "submit rejects a word that is neither an option nor KEY=VALUE" {
    run bisrde submit COUNTS=c.tsv stray
    [ "$status" -eq 2 ]
    [[ "$output" == *"expected an sbatch option or KEY=VALUE"* ]]
    [[ "$output" != *"sbatch-arg"* ]]
}

@test "submit rejects a malformed KEY" {
    run bisrde submit 1BAD=x
    [ "$status" -eq 2 ]
    [[ "$output" == *"not a KEY=VALUE setting"* ]]
}

# -------------------------------------------------- interactive and pull

@test "pull reports an image that is already present" {
    with_apptainer
    mkdir -p "$BISR_SIF_DIR"
    : > "$BISR_SIF_DIR/bisrde_$VERSION.sif"
    run bisrde pull
    [ "$status" -eq 0 ]
    [[ "$output" == *"already present: $BISR_SIF_DIR/bisrde_$VERSION.sif"* ]]
}

@test "pull creates the image through the container tool when missing" {
    with_apptainer
    unset BISR_DRY_RUN
    run bisrde pull
    [ "$status" -eq 0 ]
    [ -f "$BISR_SIF_DIR/bisrde_$VERSION.sif" ]
    [[ "$output" == *"apptainer pull $BISR_SIF_DIR/bisrde_$VERSION.sif docker://"* ]]
}

# -------------------------------------------------------------- modulefile

# Source the modulefile in tclsh with stubbed module commands and print what it
# set. $1 is the BISR_SIF_DIR value to export first, empty for none.
load_module() {
    cat > "$TMP/load.tcl" <<'EOF2'
proc prepend-path {var val} { puts "prepend-path $var $val" }
proc setenv {var val} { puts "setenv $var $val" }
proc module-whatis {text} { puts "whatis $text" }
set ModulesCurrentModulefile [lindex $argv 0]
source $ModulesCurrentModulefile
ModulesHelp
EOF2
    if [ -n "${1:-}" ]; then
        BISR_SIF_DIR="$1" /usr/bin/tclsh "$TMP/load.tcl" "$DE/modulefiles/bisrde" 2>&1
    else
        env -u BISR_SIF_DIR /usr/bin/tclsh "$TMP/load.tcl" "$DE/modulefiles/bisrde" 2>&1
    fi
}

@test "modulefile puts bin on PATH and defaults BISR_SIF_DIR to the lab directory" {
    [ -x /usr/bin/tclsh ] || skip "tclsh is not installed"
    run load_module ""
    [ "$status" -eq 0 ]
    [[ "$output" == *"prepend-path PATH $DE/bin"* ]]
    [[ "$output" == *"setenv BISR_SIF_DIR /lustre/home/mccbnfolab/containers"* ]]
    [[ "$output" == *"whatis "* ]]
    [[ "$output" == *"bisrde run"* ]]
}

@test "modulefile leaves an exported BISR_SIF_DIR alone" {
    [ -x /usr/bin/tclsh ] || skip "tclsh is not installed"
    run load_module "$TMP/mine"
    [ "$status" -eq 0 ]
    [[ "$output" == *"prepend-path PATH $DE/bin"* ]]
    [[ "$output" != *"setenv BISR_SIF_DIR"* ]]
}

@test "modulefile starts with the Modules header" {
    [ "$(head -n 1 "$DE/modulefiles/bisrde")" = "#%Module1.0" ]
}
