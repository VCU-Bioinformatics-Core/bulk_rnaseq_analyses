package main

import (
	"fmt"
	"strings"
)

// Config holds every pipeline parameter collected from the form / env.
type Config struct {
	Counts           string
	Samplesheet      string
	Outdir           string
	Runid            string
	Annotation       string
	BRS              string
	IDType           string
	VolcanoLabels    string
	FoldChange       string
	Padj             string
	GseaRank         string // stat | log2fc
	LfcShrink        string // apeglm | normal | none
	IndepFiltering   string // yes | no
	ExcludeSamples   []string
	ExcludeGroups    []string
	IncludeContrasts []string
	ExcludeContrasts []string
}

// ToArgs returns the run_analysis.sh argument vector. This MUST stay in
// parity with run_interactive.sh's flag assembly (there is a parity test
// in the project's CLAUDE/tracker notes).
func (c Config) ToArgs() []string {
	args := []string{
		"--counts", c.Counts,
		"--samplesheet", c.Samplesheet,
		"--outdir", c.Outdir,
		"--runid", c.Runid,
		"--annotation", c.Annotation,
		"--id-type", c.IDType,
		"--volcano-labels", c.VolcanoLabels,
		"--fold-change", c.FoldChange,
		"--padj", c.Padj,
		"--gsea-rank", c.GseaRank,
		"--lfc-shrink", c.LfcShrink,
	}
	// Exact "yes" to stay byte-identical with run_interactive.sh's test.
	if strings.TrimSpace(c.IndepFiltering) == "yes" {
		args = append(args, "--independent-filtering")
	}
	if strings.TrimSpace(c.BRS) != "" {
		args = append(args, "--brs-ticket", c.BRS)
	}
	if len(c.ExcludeSamples) > 0 {
		args = append(args, "--exclude-samples", strings.Join(c.ExcludeSamples, ","))
	}
	if len(c.ExcludeGroups) > 0 {
		args = append(args, "--exclude-groups", strings.Join(c.ExcludeGroups, ","))
	}
	if len(c.IncludeContrasts) > 0 {
		args = append(args, "--include-contrasts", strings.Join(c.IncludeContrasts, ","))
	}
	if len(c.ExcludeContrasts) > 0 {
		args = append(args, "--exclude-contrasts", strings.Join(c.ExcludeContrasts, ","))
	}
	return args
}

// CommandPreview is the human-readable run_analysis.sh invocation.
func (c Config) CommandPreview() string {
	return "bash run_analysis.sh " + strings.Join(c.ToArgs(), " ")
}

func dash(s string) string {
	if strings.TrimSpace(s) == "" {
		return "—"
	}
	return s
}

func dashList(s []string) string {
	if len(s) == 0 {
		return "—"
	}
	return strings.Join(s, ", ")
}

// Summary renders a config summary block.
func (c Config) Summary() string {
	rows := [][2]string{
		{"Counts", c.Counts},
		{"Samplesheet", c.Samplesheet},
		{"Annotation", c.Annotation},
		{"Run ID", c.Runid},
		{"Output dir", c.Outdir},
		{"BRS ticket", dash(c.BRS)},
		{"ID type", c.IDType},
		{"Volcano labels", c.VolcanoLabels},
		{"Fold-change", c.FoldChange},
		{"Adj p-value", c.Padj},
		{"GSEA ranking", c.GseaRank},
		{"LFC shrinkage", c.LfcShrink},
		{"Indep. filtering", c.IndepFiltering},
		{"Exclude groups", dashList(c.ExcludeGroups)},
		{"Exclude samples", dashList(c.ExcludeSamples)},
		{"Include contrasts", dashList(c.IncludeContrasts)},
		{"Exclude contrasts", dashList(c.ExcludeContrasts)},
	}
	var b strings.Builder
	b.WriteString("Run configuration\n\n")
	for _, r := range rows {
		b.WriteString(fmt.Sprintf("%s  %s\n", keyStyle.Render(fmt.Sprintf("%-18s", r[0])), r[1]))
	}
	return strings.TrimRight(b.String(), "\n")
}
