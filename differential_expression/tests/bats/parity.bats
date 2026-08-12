#!/usr/bin/env bats
#
# Launcher parity + samplesheet handling.
#
# run_interactive.sh (bash) and tui/bisrde-tui (Go) must assemble a
# BYTE-IDENTICAL run_analysis.sh argument vector — that is the maintained
# contract from ADR-0005, and the two implementations are independent, so
# they can silently drift. The Go side has tui/samplesheet_test.go; this is
# the bash side plus the cross-implementation assertions.
#
# Run with:  just test-sh     (or: bats tests/bats)

setup() {
    DE="$( cd "$( dirname "$BATS_TEST_FILENAME" )/../.." && pwd )"
    export DE
    TUI="$DE/tui/bisrde-tui"
    export TUI
    # Common non-interactive config. Every var is set so no prompt can block.
    export BISR_FRONTEND=bash BISR_ASSUME_YES=1 \
           BISR_COUNTS=assets/example_counts.tsv \
           BISR_SAMPLESHEET=assets/example_samplesheet.csv \
           BISR_ANNOTATION=mouse BISR_RUNID=parity \
           BISR_OUTDIR=./results BISR_IDTYPE=ensembl
}

# Assemble via the bash launcher.
bash_cmd() {
    ( cd "$DE" && bash run_interactive.sh --print-cmd 2>/dev/null \
        | grep '^bash run_analysis.sh' )
}

# Assemble via the Go TUI (builds on demand).
tui_cmd() {
    if [ ! -x "$TUI" ]; then
        ( cd "$DE/tui" && go build -o bisrde-tui . ) || return 1
    fi
    ( cd "$DE" && "$TUI" --print-cmd 2>/dev/null | grep '^bash run_analysis.sh' )
}

@test "bash launcher emits a run_analysis.sh command" {
    run bash_cmd
    [ "$status" -eq 0 ]
    [[ "$output" == bash\ run_analysis.sh\ * ]]
}

@test "bash and Go TUI assemble byte-identical argument vectors" {
    b="$(bash_cmd)"
    t="$(tui_cmd)"
    [ -n "$b" ]
    [ "$b" = "$t" ]
}

@test "parity holds with an optional BRS ticket set" {
    export BISR_BRS=BRS-1234
    b="$(bash_cmd)"; t="$(tui_cmd)"
    [ "$b" = "$t" ]
    [[ "$b" == *"--brs-ticket BRS-1234"* ]]
}

@test "parity holds with non-default thresholds and volcano labels" {
    export BISR_FOLD_CHANGE=2 BISR_PADJ=0.01 BISR_VOLCANO_LABELS=25
    b="$(bash_cmd)"; t="$(tui_cmd)"
    [ "$b" = "$t" ]
    [[ "$b" == *"--fold-change 2"* ]]
    [[ "$b" == *"--padj 0.01"* ]]
    [[ "$b" == *"--volcano-labels 25"* ]]
}

@test "parity holds with contrast include/exclude selection" {
    export BISR_EXCLUDE_CONTRASTS=cmtf_vs_cmtm
    b="$(bash_cmd)"; t="$(tui_cmd)"
    [ "$b" = "$t" ]
    [[ "$b" == *"--exclude-contrasts cmtf_vs_cmtm"* ]]
}

@test "every required flag is present in the assembled command" {
    b="$(bash_cmd)"
    for flag in --counts --samplesheet --outdir --runid --annotation --id-type; do
        [[ "$b" == *"$flag "* ]] || { echo "missing $flag in: $b"; return 1; }
    done
}

@test "--print-cmd does not execute the pipeline" {
    # A results dir must not be created by a dry-run assembly.
    tmp="$(mktemp -d)"
    BISR_OUTDIR="$tmp/should_not_exist" run bash_cmd
    [ "$status" -eq 0 ]
    [ ! -d "$tmp/should_not_exist" ]
    rm -rf "$tmp"
}

@test "shellcheck reports no findings in any launcher script" {
    if ! command -v shellcheck >/dev/null 2>&1; then
        skip "shellcheck not installed"
    fi
    run bash -c "cd '$DE' && shellcheck -x ./*.sh"
    [ "$status" -eq 0 ]
}

@test "samplesheet header variants are accepted (Sample_ID/Group_ID)" {
    # Mirrors the R-side test-samplesheet-columns.R: a header variant must not
    # silently produce zero contrasts.
    tmp="$(mktemp -d)"
    cat > "$tmp/ss.csv" <<'CSV'
Sample_ID,Group_ID,a_vs_b
s1,A,1
s2,B,0
CSV
    run bash -c "cd '$DE' && Rscript -e '
      suppressMessages(devtools::load_all(\"bisrDE\", quiet=TRUE))
      ss <- read_samplesheet(\"$tmp/ss.csv\")
      stopifnot(all(c(\"SampleID\",\"GroupID\") %in% colnames(ss)))
      cat(length(parse_contrasts(ss)), \"\n\")' 2>/dev/null | tail -1"
    rm -rf "$tmp"
    [ "$status" -eq 0 ]
    [[ "$output" == *"1"* ]]
}
