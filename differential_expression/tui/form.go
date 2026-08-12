package main

import (
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/charmbracelet/huh"
)

func envOr(key, def string) string {
	if v := os.Getenv(key); v != "" {
		return v
	}
	return def
}

func envList(key string) []string {
	v := os.Getenv(key)
	if strings.TrimSpace(v) == "" {
		return nil
	}
	var out []string
	for _, p := range strings.Split(v, ",") {
		if p = strings.TrimSpace(p); p != "" {
			out = append(out, p)
		}
	}
	return out
}

// defaultConfig seeds a Config from BISR_* env vars (or sensible defaults).
func defaultConfig() Config {
	return Config{
		Counts:           envOr("BISR_COUNTS", "assets/example_counts.tsv"),
		Samplesheet:      envOr("BISR_SAMPLESHEET", "assets/example_samplesheet.csv"),
		Outdir:           envOr("BISR_OUTDIR", "./results"),
		Runid:            envOr("BISR_RUNID", "run_"+time.Now().Format("20060102_150405")),
		Annotation:       envOr("BISR_ANNOTATION", "mouse"),
		BRS:              envOr("BISR_BRS", ""),
		IDType:           envOr("BISR_IDTYPE", "ensembl"),
		VolcanoLabels:    envOr("BISR_VOLCANO_LABELS", "10"),
		FoldChange:       envOr("BISR_FOLD_CHANGE", "1.5"),
		Padj:             envOr("BISR_PADJ", "0.05"),
		ExcludeSamples:   envList("BISR_EXCLUDE_SAMPLES"),
		ExcludeGroups:    envList("BISR_EXCLUDE_GROUPS"),
		IncludeContrasts: envList("BISR_INCLUDE_CONTRASTS"),
		ExcludeContrasts: envList("BISR_EXCLUDE_CONTRASTS"),
	}
}

// expandPath resolves a leading ~ and makes the path absolute.
func expandPath(p string) (string, error) {
	p = strings.TrimSpace(p)
	if strings.HasPrefix(p, "~") {
		if home, err := os.UserHomeDir(); err == nil {
			p = filepath.Join(home, strings.TrimPrefix(p, "~"))
		}
	}
	return filepath.Abs(p)
}

// pickOne offers the discovered files for one input, always with a manual
// escape hatch. Falls back to a plain text prompt when nothing was found.
func pickOne(title, desc, base string, k fileKind, dest *string) error {
	found := findCandidates(base, k, 25)
	if len(found) == 0 {
		return huh.NewForm(huh.NewGroup(
			huh.NewInput().
				Title(title).
				Description("Nothing found in " + base + " — enter a path").
				Value(dest),
		)).WithTheme(huh.ThemeCharm()).Run()
	}

	opts := make([]huh.Option[string], 0, len(found)+1)
	for _, f := range found {
		opts = append(opts, huh.NewOption(relLabel(base, f), f))
	}
	opts = append(opts, huh.NewOption(manualEntry, manualEntry))

	choice := found[0]
	if err := huh.NewForm(huh.NewGroup(
		huh.NewSelect[string]().Title(title).Description(desc).
			Options(opts...).Value(&choice),
	)).WithTheme(huh.ThemeCharm()).Run(); err != nil {
		return err
	}

	if choice == manualEntry {
		return huh.NewForm(huh.NewGroup(
			huh.NewInput().Title(title).Value(dest),
		)).WithTheme(huh.ThemeCharm()).Run()
	}
	*dest = choice
	return nil
}

// pickInputs prompts for the counts matrix and samplesheet.
func pickInputs(c *Config, base string) error {
	if err := pickOne("Counts matrix", "gene counts from nf-core/rnaseq (.tsv)",
		base, kindCounts, &c.Counts); err != nil {
		return err
	}
	return pickOne("Samplesheet", "SampleID, GroupID, then one column per contrast (.csv)",
		base, kindSamplesheet, &c.Samplesheet)
}

// BuildConfig collects parameters. In nonInteractive mode it returns the
// env/default config without prompting (used for --print-cmd and non-TTY).
func BuildConfig(nonInteractive bool) (Config, error) {
	c := defaultConfig()
	if nonInteractive {
		return c, nil
	}

	// Stage 0 — project directory, then offer the input files found inside it
	// so nobody has to type a full path from memory.
	base := envOr("BISR_BASE_DIR", "")
	if base == "" {
		base, _ = os.Getwd()
	}
	baseForm := huh.NewForm(huh.NewGroup(
		huh.NewInput().
			Title("Project directory").
			Description("Searched for counts / samplesheet files").
			Value(&base),
	)).WithTheme(huh.ThemeCharm())
	if err := baseForm.Run(); err != nil {
		return c, err
	}
	if expanded, err := expandPath(base); err == nil {
		base = expanded
	}
	if err := pickInputs(&c, base); err != nil {
		return c, err
	}

	// Stage 1 — core parameters. Counts/samplesheet are already chosen above,
	// but stay editable here.
	core := huh.NewForm(
		huh.NewGroup(
			huh.NewInput().Title("Counts TSV").Value(&c.Counts),
			huh.NewInput().Title("Samplesheet CSV").Value(&c.Samplesheet),
			huh.NewSelect[string]().Title("Annotation").
				Options(huh.NewOptions("mouse", "human")...).
				Value(&c.Annotation),
			huh.NewInput().Title("Run ID").Value(&c.Runid),
			huh.NewInput().Title("Output directory").Value(&c.Outdir),
			huh.NewInput().Title("BRS ticket (optional)").Value(&c.BRS),
			huh.NewSelect[string]().Title("Gene ID type").
				Options(huh.NewOptions("ensembl", "entrez", "symbol")...).
				Value(&c.IDType),
			huh.NewInput().Title("Volcano labels (top N genes)").Value(&c.VolcanoLabels),
			huh.NewInput().Title("Fold-change cutoff (e.g. 2 = 2-fold)").Value(&c.FoldChange),
			huh.NewInput().Title("Adjusted p-value (FDR) cutoff").Value(&c.Padj),
		),
	).WithTheme(huh.ThemeCharm())
	if err := core.Run(); err != nil {
		return c, err
	}

	// Stage 2 — selection menus, derived from the chosen samplesheet.
	ss, err := ReadSamplesheet(c.Samplesheet)
	if err != nil || len(ss.Header) == 0 {
		return c, nil // samplesheet unreadable -> no filtering offered
	}
	groups := ss.Groups()
	samples := ss.Samples()
	contrasts := ss.Contrasts()

	contrastMode := "all"
	sel := huh.NewForm(
		huh.NewGroup(
			huh.NewMultiSelect[string]().Title("Exclude groups (optional)").
				Options(huh.NewOptions(groups...)...).Value(&c.ExcludeGroups),
			huh.NewMultiSelect[string]().Title("Exclude samples (optional)").
				Options(huh.NewOptions(samples...)...).Value(&c.ExcludeSamples),
			huh.NewSelect[string]().Title("Contrasts to run").
				Options(
					huh.NewOption("all", "all"),
					huh.NewOption("include only a subset", "include"),
					huh.NewOption("exclude a subset", "exclude"),
				).Value(&contrastMode),
		),
	).WithTheme(huh.ThemeCharm())
	if err := sel.Run(); err != nil {
		return c, err
	}

	switch contrastMode {
	case "include":
		f := huh.NewForm(huh.NewGroup(
			huh.NewMultiSelect[string]().Title("Include these contrasts").
				Options(huh.NewOptions(contrasts...)...).Value(&c.IncludeContrasts),
		)).WithTheme(huh.ThemeCharm())
		if err := f.Run(); err != nil {
			return c, err
		}
	case "exclude":
		f := huh.NewForm(huh.NewGroup(
			huh.NewMultiSelect[string]().Title("Exclude these contrasts").
				Options(huh.NewOptions(contrasts...)...).Value(&c.ExcludeContrasts),
		)).WithTheme(huh.ThemeCharm())
		if err := f.Run(); err != nil {
			return c, err
		}
	}

	// Stage 4 — sample display names (#2b). When the sheet has no usable
	// DisplayName column, offer per-group labels (per-sample doesn't scale);
	// otherwise the R side auto-derives "<Group> <n>" with a warning. Skipped
	// entirely in non-interactive mode (early return above), so --print-cmd
	// parity with run_interactive.sh is preserved.
	if len(groups) > 0 && !ss.HasDisplayNames() {
		dnMode := "auto"
		dn := huh.NewForm(huh.NewGroup(
			huh.NewSelect[string]().
				Title("Sample display names").
				Description("No DisplayName column found. Tip: add one to the samplesheet for full control.").
				Options(
					huh.NewOption("Auto-derive from group + replicate (recommended)", "auto"),
					huh.NewOption("Enter a label per group", "group"),
				).Value(&dnMode),
		)).WithTheme(huh.ThemeCharm())
		if err := dn.Run(); err != nil {
			return c, err
		}
		if dnMode == "group" {
			labelVals := make([]string, len(groups)) // stable backing array for the pointers
			fields := make([]huh.Field, len(groups))
			for i, g := range groups {
				labelVals[i] = g // default label = the group name
				fields[i] = huh.NewInput().Title("Label for group '" + g + "'").Value(&labelVals[i])
			}
			lf := huh.NewForm(huh.NewGroup(fields...)).WithTheme(huh.ThemeCharm())
			if err := lf.Run(); err != nil {
				return c, err
			}
			labels := make(map[string]string, len(groups))
			for i, g := range groups {
				labels[g] = labelVals[i]
			}
			if path, err := ss.WriteDisplaySamplesheet(labels); err == nil {
				c.Samplesheet = path
			}
		}
	}
	return c, nil
}
