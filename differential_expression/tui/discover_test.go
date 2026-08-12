package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// buildTree creates files (relative paths) under a fresh temp dir.
func buildTree(t *testing.T, paths ...string) string {
	t.Helper()
	root := t.TempDir()
	for _, p := range paths {
		full := filepath.Join(root, p)
		if err := os.MkdirAll(filepath.Dir(full), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(full, []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	return root
}

func names(base string, got []string) []string {
	out := make([]string, 0, len(got))
	for _, g := range got {
		out = append(out, relLabel(base, g))
	}
	return out
}

func TestFindCountsPrefersLikelyNames(t *testing.T) {
	root := buildTree(t,
		"notes.txt",
		"salmon.merged.gene_counts.tsv",
		"random.tsv",
	)
	got := names(root, findCandidates(root, kindCounts, 25))
	if len(got) == 0 {
		t.Fatal("no counts candidates found")
	}
	if got[0] != "salmon.merged.gene_counts.tsv" {
		t.Errorf("best guess = %q, want the merged gene counts file; all=%v", got[0], got)
	}
}

func TestFindSamplesheetPrefersLikelyNames(t *testing.T) {
	root := buildTree(t, "misc.csv", "deg_human_ss.csv")
	got := names(root, findCandidates(root, kindSamplesheet, 25))
	if len(got) == 0 {
		t.Fatal("no samplesheet candidates found")
	}
	if got[0] != "deg_human_ss.csv" {
		t.Errorf("best guess = %q, want deg_human_ss.csv; all=%v", got[0], got)
	}
}

// The important one: a previous run's output must not pollute the pickers.
func TestPipelineOutputIsExcluded(t *testing.T) {
	root := buildTree(t,
		"deg_ss.csv",
		"salmon.merged.gene_counts.tsv",
		"results/data/de_data/DESeq2_A_vs_B.csv",
		"results/data/de_data/normalizedCounts_tmm2026-01-01.csv",
		"results/figures/qc/plot.csv",
		"results/logs/20260101_session.csv",
		"renv/library/pkg/data.csv",
		".git/objects/thing.csv",
	)

	got := names(root, findCandidates(root, kindSamplesheet, 50))
	for _, g := range got {
		for _, bad := range []string{"de_data", "figures", "logs", "renv", ".git"} {
			if strings.Contains(g, bad) {
				t.Errorf("output/noise file leaked into candidates: %q", g)
			}
		}
	}
	if len(got) != 1 || got[0] != "deg_ss.csv" {
		t.Errorf("candidates = %v, want just deg_ss.csv", got)
	}
}

func TestExtensionFiltering(t *testing.T) {
	root := buildTree(t, "a.tsv", "b.txt", "c.csv", "d.png", "e.rds")

	counts := names(root, findCandidates(root, kindCounts, 25))
	for _, c := range counts {
		if !strings.HasSuffix(c, ".tsv") && !strings.HasSuffix(c, ".txt") {
			t.Errorf("counts candidate with wrong extension: %q", c)
		}
	}
	sheets := names(root, findCandidates(root, kindSamplesheet, 25))
	for _, s := range sheets {
		if !strings.HasSuffix(s, ".csv") {
			t.Errorf("samplesheet candidate with wrong extension: %q", s)
		}
	}
	if len(counts) != 2 || len(sheets) != 1 {
		t.Errorf("counts=%v sheets=%v", counts, sheets)
	}
}

func TestScanDepthIsBounded(t *testing.T) {
	root := buildTree(t,
		"shallow.csv",
		"a/b/deep.csv",
		"a/b/c/d/e/too_deep.csv",
	)
	got := names(root, findCandidates(root, kindSamplesheet, 50))
	for _, g := range got {
		if strings.Contains(g, "too_deep") {
			t.Errorf("walked past the depth bound: %v", got)
		}
	}
	if len(got) == 0 {
		t.Fatal("depth bound excluded everything")
	}
}

func TestShallowerPathsRankFirstOnTies(t *testing.T) {
	root := buildTree(t, "sub/dir/ss.csv", "ss.csv")
	got := names(root, findCandidates(root, kindSamplesheet, 25))
	if got[0] != "ss.csv" {
		t.Errorf("expected the shallower file first, got %v", got)
	}
}

func TestFindCandidatesHandlesMissingDir(t *testing.T) {
	if got := findCandidates(filepath.Join(t.TempDir(), "nope"), kindCounts, 10); got != nil {
		t.Errorf("expected nil for a missing directory, got %v", got)
	}
	// A file (not a directory) must also be handled gracefully.
	f := filepath.Join(t.TempDir(), "f.tsv")
	os.WriteFile(f, []byte("x"), 0o644)
	if got := findCandidates(f, kindCounts, 10); got != nil {
		t.Errorf("expected nil when base is a file, got %v", got)
	}
}

func TestLimitIsRespected(t *testing.T) {
	paths := make([]string, 0, 30)
	for i := 0; i < 30; i++ {
		paths = append(paths, string(rune('a'+i%26))+string(rune('a'+i/26))+".csv")
	}
	root := buildTree(t, paths...)
	if got := findCandidates(root, kindSamplesheet, 5); len(got) != 5 {
		t.Errorf("limit ignored: got %d candidates", len(got))
	}
}

func TestRelLabelFallsBackToAbsolute(t *testing.T) {
	if got := relLabel("/a/b", "/a/b/c.csv"); got != "c.csv" {
		t.Errorf("relLabel = %q, want c.csv", got)
	}
	// Outside the base: keep the full path rather than a ../.. soup.
	if got := relLabel("/a/b", "/x/y/z.csv"); got != "/x/y/z.csv" {
		t.Errorf("relLabel = %q, want the absolute path", got)
	}
}
