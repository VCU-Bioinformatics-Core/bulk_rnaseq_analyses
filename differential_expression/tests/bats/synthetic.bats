#!/usr/bin/env bats
#
# The seeded synthetic counts generator (tests/synthetic_counts.R) and the
# shrunken-column comparison in compare_runs.R.
#
# The generator needs AnnotationDbi and org.Hs.eg.db, so these tests call the
# Rscript on PATH and inherit the environment: the renv library (export
# RENV_PATHS_LIBRARY when running from a git worktree) or the conda
# environment must be available.
#
# Run with:  just test-sh     (or: bats tests/bats)

setup() {
    DE="$( cd "$( dirname "$BATS_TEST_FILENAME" )/../.." && pwd -P )"
    export DE
    TMP="$( cd "$(mktemp -d)" && pwd -P )"
    export TMP
}

teardown() {
    rm -rf "$TMP"
}

# Two fabricated run directories with one comparison, 12 genes. Run A keeps
# the shrunken columns as given; run B gets the same table.
fabricate_run() {
    local dir="$1" shrunken="${2:-yes}"
    mkdir -p "$dir/data/de_data"
    local f="$dir/data/de_data/DESeq2_x.csv"
    if [ "$shrunken" = yes ]; then
        echo "ENSEMBL_ID,log2FoldChange,padj,log2FC_shrunken,lfcSE_shrunken" > "$f"
    else
        echo "ENSEMBL_ID,log2FoldChange,padj" > "$f"
    fi
    local i
    for i in $(seq 1 12); do
        local lfc padj shr se
        lfc=$(awk -v i="$i" 'BEGIN { printf "%.4f", (i - 6.5) * 0.7 }')
        padj=$(awk -v i="$i" 'BEGIN { printf "%.6f", i * 0.004 }')
        shr=$(awk -v i="$i" 'BEGIN { printf "%.4f", (i - 6.5) * 0.5 }')
        se=$(awk -v i="$i" 'BEGIN { printf "%.4f", 0.1 + i / 100 }')
        if [ "$shrunken" = yes ]; then
            echo "G$i,$lfc,$padj,$shr,$se" >> "$f"
        else
            echo "G$i,$lfc,$padj" >> "$f"
        fi
    done
}

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}

@test "generator writes the four files with the requested sizes" {
    run Rscript "$DE/tests/synthetic_counts.R" --out "$TMP/d" --genes 300 --per-group 2
    [ "$status" -eq 0 ]
    for f in counts.tsv samplesheet.csv truth.csv manifest.txt; do
        [ -f "$TMP/d/$f" ]
    done
    # Header plus 300 genes; gene_id plus 4 samples.
    [ "$(wc -l < "$TMP/d/counts.tsv" | tr -d ' ')" -eq 301 ]
    [ "$(head -1 "$TMP/d/counts.tsv" | awk -F'\t' '{print NF}')" -eq 5 ]
    [ "$(awk -F, 'NR > 1 && $2 == "TRUE"' "$TMP/d/truth.csv" | wc -l | tr -d ' ')" -eq 30 ]
    [ "$(head -1 "$TMP/d/samplesheet.csv")" = "SampleID,GroupID,treated_vs_control" ]
    grep -q '^seed: 20261015$' "$TMP/d/manifest.txt"
    grep -q '^org.Hs.eg.db: ' "$TMP/d/manifest.txt"
}

@test "generator is deterministic for one seed and differs across seeds" {
    Rscript "$DE/tests/synthetic_counts.R" --out "$TMP/a" --genes 300 --per-group 2
    Rscript "$DE/tests/synthetic_counts.R" --out "$TMP/b" --genes 300 --per-group 2
    Rscript "$DE/tests/synthetic_counts.R" --out "$TMP/c" --genes 300 --per-group 2 --seed 7
    [ "$(md5_of "$TMP/a/counts.tsv")" = "$(md5_of "$TMP/b/counts.tsv")" ]
    [ "$(md5_of "$TMP/a/truth.csv")" = "$(md5_of "$TMP/b/truth.csv")" ]
    [ "$(md5_of "$TMP/a/counts.tsv")" != "$(md5_of "$TMP/c/counts.tsv")" ]
}

@test "generator requires --out" {
    run Rscript "$DE/tests/synthetic_counts.R" --genes 10
    [ "$status" -eq 3 ]
}

@test "compare_runs prints the shrunken lines and passes identical runs" {
    fabricate_run "$TMP/A"
    fabricate_run "$TMP/B"
    run Rscript "$DE/compare_runs.R" --a "$TMP/A" --b "$TMP/B"
    [ "$status" -eq 0 ]
    [[ "$output" == *"- log2FC_shrunken: Spearman=1.0000 Pearson=1.0000 median|Δ|=0.0000 max|Δ|=0.0000"* ]]
    [[ "$output" == *"- lfcSE_shrunken: max|Δ|=0.0000"* ]]
    [[ "$output" == *"## Overall: PASS"* ]]
}

@test "compare_runs fails when only the shrunken column differs" {
    fabricate_run "$TMP/A"
    fabricate_run "$TMP/B"
    # Reverse the shrunken estimates; log2FoldChange and padj stay identical.
    f="$TMP/B/data/de_data/DESeq2_x.csv"
    awk -F, -v OFS=, 'NR == FNR { if (FNR > 1) v[FNR] = $4; n = FNR; next }
                      FNR == 1 { print; next }
                      { $4 = v[n + 2 - FNR]; print }' "$f" "$f" > "$f.new"
    mv "$f.new" "$f"
    run Rscript "$DE/compare_runs.R" --a "$TMP/A" --b "$TMP/B"
    [ "$status" -eq 2 ]
    [[ "$output" == *"- log2FC: Spearman=1.0000"* ]]
    [[ "$output" == *"- log2FC_shrunken: Spearman=-1.0000"* ]]
    [[ "$output" == *"## Overall: FAIL"* ]]
}

@test "compare_runs notes a one-sided shrunken column and keeps the verdict" {
    fabricate_run "$TMP/A"
    fabricate_run "$TMP/B" no
    run Rscript "$DE/compare_runs.R" --a "$TMP/A" --b "$TMP/B"
    [ "$status" -eq 0 ]
    [[ "$output" == *"- shrunken columns: only in A, not compared"* ]]
    [[ "$output" != *"- log2FC_shrunken:"* ]]
    [[ "$output" == *"## Overall: PASS"* ]]
}
