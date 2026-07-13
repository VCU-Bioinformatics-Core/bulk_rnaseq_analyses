package main

import (
	"os"
	"testing"
)

func writeTempCSV(t *testing.T, content string) string {
	t.Helper()
	f, err := os.CreateTemp("", "ss_*.csv")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := f.WriteString(content); err != nil {
		t.Fatal(err)
	}
	f.Close()
	return f.Name()
}

func TestHasDisplayNames(t *testing.T) {
	cases := []struct {
		name    string
		content string
		want    bool
	}{
		{"no column", "SampleID,GroupID,c1\nS1,A,1\nS2,B,0\n", false},
		{"all blank", "SampleID,GroupID,DisplayName,c1\nS1,A,,1\nS2,B,,0\n", false},
		{"has value", "SampleID,GroupID,DisplayName,c1\nS1,A,Ctrl 1,1\nS2,B,,0\n", true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			p := writeTempCSV(t, tc.content)
			defer os.Remove(p)
			ss, err := ReadSamplesheet(p)
			if err != nil {
				t.Fatal(err)
			}
			if got := ss.HasDisplayNames(); got != tc.want {
				t.Errorf("HasDisplayNames() = %v, want %v", got, tc.want)
			}
		})
	}
}

func TestWriteDisplaySamplesheet(t *testing.T) {
	src := writeTempCSV(t, "SampleID,GroupID,c1\nS1,pos,1\nS2,pos,1\nS3,neg,0\nS4,pos,1\n")
	defer os.Remove(src)
	ss, err := ReadSamplesheet(src)
	if err != nil {
		t.Fatal(err)
	}

	out, err := ss.WriteDisplaySamplesheet(map[string]string{"pos": "HPV+", "neg": "HPV-"})
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(out)

	got, err := ReadSamplesheet(out)
	if err != nil {
		t.Fatal(err)
	}
	di := got.colIndex("DisplayName")
	if di < 0 {
		t.Fatal("DisplayName column not written")
	}
	// Replicates numbered within each group, in row order.
	want := []string{"HPV+ 1", "HPV+ 2", "HPV- 1", "HPV+ 3"}
	for i, row := range got.Rows {
		if row[di] != want[i] {
			t.Errorf("row %d DisplayName = %q, want %q", i, row[di], want[i])
		}
	}
}

func TestWriteDisplaySamplesheet_FallbackToGroup(t *testing.T) {
	src := writeTempCSV(t, "SampleID,GroupID,c1\nS1,pos,1\nS2,neg,0\n")
	defer os.Remove(src)
	ss, _ := ReadSamplesheet(src)
	out, err := ss.WriteDisplaySamplesheet(map[string]string{}) // no labels -> GroupID
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(out)
	got, _ := ReadSamplesheet(out)
	di := got.colIndex("DisplayName")
	if got.Rows[0][di] != "pos 1" || got.Rows[1][di] != "neg 1" {
		t.Errorf("fallback labels = %q, %q; want \"pos 1\", \"neg 1\"",
			got.Rows[0][di], got.Rows[1][di])
	}
}
