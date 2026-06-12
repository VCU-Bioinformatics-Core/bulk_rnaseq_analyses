package main

import (
	"os"
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
		ExcludeSamples:   envList("BISR_EXCLUDE_SAMPLES"),
		ExcludeGroups:    envList("BISR_EXCLUDE_GROUPS"),
		IncludeContrasts: envList("BISR_INCLUDE_CONTRASTS"),
		ExcludeContrasts: envList("BISR_EXCLUDE_CONTRASTS"),
	}
}

// BuildConfig collects parameters. In nonInteractive mode it returns the
// env/default config without prompting (used for --print-cmd and non-TTY).
func BuildConfig(nonInteractive bool) (Config, error) {
	c := defaultConfig()
	if nonInteractive {
		return c, nil
	}

	// Stage 1 — core parameters.
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
	return c, nil
}
