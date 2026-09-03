// Command bisrde-tui is the Go TUI front-end (charmbracelet Phase B) for the
// bisrDE bulk RNA-seq DE pipeline. It collects parameters through a styled
// huh form (with interactive group / sample / contrast selection derived from
// the samplesheet) and hands off to run_analysis.sh — exactly the same flag
// vector that run_interactive.sh assembles.
//
// Usage:
//
//	bisrde-tui                 # full interactive TUI
//	bisrde-tui --print-cmd     # assemble + print the run_analysis.sh command
//	bisrde-tui --help
//
// Project location: defaults to the current working directory; override with
// BISR_PROJECT_DIR. The bash chooser (run_interactive.sh) sets this for you.
package main

import (
	"fmt"
	"os"

	"github.com/charmbracelet/huh"
)

func printHelp() {
	fmt.Print(`bisrde-tui — Go TUI launcher for the bisrDE pipeline (charmbracelet Phase B)

Run with no arguments for the guided TUI. It collects the counts file,
samplesheet, annotation, run ID, output dir, optional BRS ticket and ID type,
then offers interactive group / sample / contrast selection before launching.

Flags:
    --print-cmd   Assemble and print the run_analysis.sh command (no run).
    -h, --help    This help.

Env overrides (skip the matching prompt): BISR_COUNTS, BISR_SAMPLESHEET,
BISR_ANNOTATION, BISR_RUNID, BISR_OUTDIR, BISR_BRS, BISR_IDTYPE, BISR_GSEA_RANK, BISR_LFC_SHRINK, BISR_INDEPENDENT_FILTERING,
BISR_EXCLUDE_SAMPLES, BISR_EXCLUDE_GROUPS, BISR_INCLUDE_CONTRASTS,
BISR_EXCLUDE_CONTRASTS, BISR_PROJECT_DIR.
`)
}

func main() {
	printCmd := false
	for _, a := range os.Args[1:] {
		switch a {
		case "--print-cmd":
			printCmd = true
		case "-h", "--help":
			printHelp()
			return
		}
	}

	// Interactive only with a real terminal on stdin and not just printing.
	fi, _ := os.Stdin.Stat()
	isTTY := fi != nil && (fi.Mode()&os.ModeCharDevice) != 0
	nonInteractive := printCmd || !isTTY

	projectDir := envOr("BISR_PROJECT_DIR", "")
	if projectDir == "" {
		projectDir, _ = os.Getwd()
	}

	if !nonInteractive {
		fmt.Println(bannerStyle.Render(
			"Bulk RNA-Seq Differential Expression\nbisrDE · v1.6.4 · Go TUI"))
	}

	cfg, err := BuildConfig(nonInteractive)
	if err != nil {
		fmt.Fprintln(os.Stderr, errStyle.Render("aborted: "+err.Error()))
		os.Exit(1)
	}

	if printCmd {
		fmt.Println(cfg.CommandPreview())
		return
	}

	fmt.Println(paneStyle.Render(cfg.Summary()))

	if isTTY {
		ok := true
		if err := huh.NewConfirm().
			Title("Launch the pipeline with these settings?").
			Value(&ok).WithTheme(huh.ThemeCharm()).Run(); err != nil {
			fmt.Fprintln(os.Stderr, errStyle.Render("aborted: "+err.Error()))
			os.Exit(1)
		}
		if !ok {
			fmt.Println(mutedStyle.Render("Aborted — nothing was run."))
			return
		}
	}

	if err := RunPipeline(projectDir, cfg); err != nil {
		fmt.Println(errStyle.Render("✗ Pipeline failed: " + err.Error()))
		os.Exit(1)
	}
	fmt.Println(okStyle.Render("✓ Pipeline complete — output in " + cfg.Outdir))
}
