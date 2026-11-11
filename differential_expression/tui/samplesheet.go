package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"strings"
)

// metaCols are the reserved, non-contrast samplesheet columns — kept in sync
// with bisrDE's .samplesheet_meta_cols (io.R).
var metaCols = map[string]bool{
	"SampleID": true, "GroupID": true, "Exclude": true, "DisplayName": true,
}

// Samplesheet is a parsed CSV (header + rows).
type Samplesheet struct {
	Header []string
	Rows   [][]string
}

// ReadSamplesheet parses a comma-delimited samplesheet.
func ReadSamplesheet(path string) (*Samplesheet, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.FieldsPerRecord = -1 // tolerate ragged rows
	recs, err := r.ReadAll()
	if err != nil {
		return nil, err
	}
	if len(recs) == 0 {
		return &Samplesheet{}, nil
	}
	return &Samplesheet{Header: recs[0], Rows: recs[1:]}, nil
}

func (s *Samplesheet) colIndex(name string) int {
	for i, h := range s.Header {
		if strings.TrimSpace(h) == name {
			return i
		}
	}
	return -1
}

// Samples returns the SampleID column (default: first column).
func (s *Samplesheet) Samples() []string {
	idx := s.colIndex("SampleID")
	if idx < 0 {
		idx = 0
	}
	return s.column(idx, false)
}

// Groups returns the unique GroupID values (default: second column).
func (s *Samplesheet) Groups() []string {
	idx := s.colIndex("GroupID")
	if idx < 0 {
		idx = 1
	}
	return s.column(idx, true)
}

// Contrasts returns the header columns that are not reserved metadata.
func (s *Samplesheet) Contrasts() []string {
	var out []string
	for _, h := range s.Header {
		ht := strings.TrimSpace(h)
		if ht != "" && !metaCols[ht] {
			out = append(out, ht)
		}
	}
	return out
}

func (s *Samplesheet) column(idx int, unique bool) []string {
	seen := map[string]bool{}
	var out []string
	for _, row := range s.Rows {
		if idx < 0 || idx >= len(row) {
			continue
		}
		v := strings.TrimSpace(row[idx])
		if v == "" {
			continue
		}
		if unique {
			if seen[v] {
				continue
			}
			seen[v] = true
		}
		out = append(out, v)
	}
	return out
}

// HasDisplayNames reports whether the sheet carries a DisplayName column with
// at least one non-blank value. Mirrors ss_has_display in run_interactive.sh;
// drives the #2b per-group prompt.
func (s *Samplesheet) HasDisplayNames() bool {
	idx := s.colIndex("DisplayName")
	if idx < 0 {
		return false
	}
	for _, row := range s.Rows {
		if idx < len(row) && strings.TrimSpace(row[idx]) != "" {
			return true
		}
	}
	return false
}

// WriteDisplaySamplesheet writes a working copy of the samplesheet with a
// DisplayName column filled from a group->label map, numbering replicates
// within each group ("<label> <n>", in row order). A blank/absent label falls
// back to the GroupID. Mirrors write_display_samplesheet in run_interactive.sh.
// Returns the path to the working copy (a temp file).
func (s *Samplesheet) WriteDisplaySamplesheet(labels map[string]string) (string, error) {
	gIdx := s.colIndex("GroupID")
	if gIdx < 0 {
		gIdx = 1
	}
	dIdx := s.colIndex("DisplayName") // -1 -> append a new column

	header := append([]string{}, s.Header...)
	if dIdx < 0 {
		header = append(header, "DisplayName")
	}

	out := [][]string{header}
	counts := map[string]int{}
	for _, row := range s.Rows {
		r := append([]string{}, row...)
		group := ""
		if gIdx >= 0 && gIdx < len(r) {
			group = strings.TrimSpace(r[gIdx])
		}
		counts[group]++
		label := group
		if l, ok := labels[group]; ok && strings.TrimSpace(l) != "" {
			label = strings.TrimSpace(l)
		}
		dn := fmt.Sprintf("%s %d", label, counts[group])
		if dIdx < 0 {
			r = append(r, dn)
		} else {
			for len(r) <= dIdx {
				r = append(r, "")
			}
			r[dIdx] = dn
		}
		out = append(out, r)
	}

	f, err := os.CreateTemp("", "bisr_ss_*.csv")
	if err != nil {
		return "", err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	if err := w.WriteAll(out); err != nil {
		return "", err
	}
	if err := w.Error(); err != nil {
		return "", err
	}
	return f.Name(), nil
}
