package main

import (
	"encoding/csv"
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
