package main

import (
	"fmt"
	"os"
	"os/exec"
)

// RunPipeline shells out to run_analysis.sh and streams its output natively.
//
// We intentionally do NOT capture the R pipeline's stdout into a bubbletea
// viewport: the pipeline emits rich `cli` progress (banners, progress bars,
// colour) that is best shown live and unmodified — the same reasoning as
// run_interactive.sh. The TUI's value is the styled huh form + lipgloss
// panels around the run, not re-rendering the pipeline's own output. A live
// multi-pane log viewport is a possible future enhancement (see tui/README).
func RunPipeline(projectDir string, c Config) error {
	fmt.Println(mutedStyle.Render("Launching… (full pipeline output follows)"))
	fmt.Println(mutedStyle.Render(c.CommandPreview()))
	fmt.Println()

	cmd := exec.Command("bash", append([]string{"run_analysis.sh"}, c.ToArgs()...)...)
	cmd.Dir = projectDir
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	cmd.Stdin = os.Stdin
	return cmd.Run()
}
