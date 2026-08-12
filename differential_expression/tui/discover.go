package main

import (
	"io/fs"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// Input discovery.
//
// Typing a full path to a counts matrix is the least pleasant part of starting
// a run. Instead we ask for a project directory once and offer the plausible
// files found inside it.

type fileKind int

const (
	kindCounts fileKind = iota
	kindSamplesheet
)

// manualEntry is the sentinel option offered alongside discovered files so a
// path outside the scanned tree is always reachable.
const manualEntry = "→ enter a path manually…"

// maxScanDepth bounds the walk. Inputs generally sit at the project root or one
// or two levels down (e.g. human/run2/salmon.merged.gene_counts.tsv); walking
// deeper mostly finds pipeline output.
const maxScanDepth = 3

// skipDirs are never descended into. The pipeline's OWN output matters most
// here: without this a previous run's data/de_data/DESeq2_*.csv files would
// swamp the samplesheet picker, and figures/logs are pure noise.
var skipDirs = map[string]bool{
	"de_data": true, "gsea_data": true, "kegg_data": true,
	"reactome_data": true, "hallmark_data": true,
	"figures": true, "logs": true, "renv": true, "work": true,
	".git": true, ".quarto": true, "node_modules": true,
	".Rproj.user": true, "__pycache__": true, ".venv": true,
}

// candidateExts returns the extensions worth offering for a kind.
func candidateExts(k fileKind) []string {
	if k == kindCounts {
		return []string{".tsv", ".txt"}
	}
	return []string{".csv"}
}

// score ranks a filename by how likely it is to be the file we're asking for.
// Higher is better; ties fall back to path depth then alphabetical order.
func score(path string, k fileKind) int {
	name := strings.ToLower(filepath.Base(path))
	s := 0
	if k == kindCounts {
		for _, hint := range []string{"counts", "merged", "gene", "matrix"} {
			if strings.Contains(name, hint) {
				s += 10
			}
		}
		if strings.HasSuffix(name, ".tsv") {
			s += 5
		}
	} else {
		for _, hint := range []string{"samplesheet", "sample_sheet", "sheet", "samples", "deg", "_ss", "ss."} {
			if strings.Contains(name, hint) {
				s += 10
			}
		}
		// A previous run's DE output is never the samplesheet.
		if strings.HasPrefix(name, "deseq2_") || strings.Contains(name, "normalizedcounts") {
			s -= 50
		}
	}
	return s
}

// findCandidates walks base (bounded depth, skipping output directories) and
// returns plausible input files, best guess first, capped at limit.
func findCandidates(base string, k fileKind, limit int) []string {
	info, err := os.Stat(base)
	if err != nil || !info.IsDir() {
		return nil
	}
	exts := candidateExts(k)
	baseDepth := strings.Count(filepath.Clean(base), string(os.PathSeparator))

	var found []string
	_ = filepath.WalkDir(base, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return nil // unreadable entries are skipped, never fatal
		}
		if d.IsDir() {
			if p == base {
				return nil
			}
			name := d.Name()
			if skipDirs[name] || strings.HasPrefix(name, ".") {
				return filepath.SkipDir
			}
			if strings.Count(filepath.Clean(p), string(os.PathSeparator))-baseDepth >= maxScanDepth {
				return filepath.SkipDir
			}
			return nil
		}
		ext := strings.ToLower(filepath.Ext(p))
		for _, want := range exts {
			if ext == want {
				found = append(found, p)
				break
			}
		}
		return nil
	})

	sort.SliceStable(found, func(i, j int) bool {
		si, sj := score(found[i], k), score(found[j], k)
		if si != sj {
			return si > sj
		}
		di := strings.Count(found[i], string(os.PathSeparator))
		dj := strings.Count(found[j], string(os.PathSeparator))
		if di != dj {
			return di < dj // shallower first
		}
		return found[i] < found[j]
	})

	if limit > 0 && len(found) > limit {
		found = found[:limit]
	}
	return found
}

// relLabel shortens a path for display, relative to base when possible.
func relLabel(base, path string) string {
	if rel, err := filepath.Rel(base, path); err == nil && !strings.HasPrefix(rel, "..") {
		return rel
	}
	return path
}
