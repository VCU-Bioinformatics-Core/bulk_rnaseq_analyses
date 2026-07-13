package main

import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"os/signal"
	"syscall"

	"github.com/creack/pty"
	"github.com/mattn/go-isatty"
)

// RunPipeline shells out to run_analysis.sh and streams its output natively.
//
// When stdout is a real terminal we run the child under a pseudo-terminal
// (creack/pty) so R's `cli` sees a TTY and renders live progress bars.
// Otherwise cli autodetects a non-TTY stream and falls back to one-shot
// milestone prints (the "67% -> 100%, printed once" symptom). Note cli emits
// progress on *stderr* when the pipeline's session-log sink is active, and
// pty.Start wires the child's stdin, stdout AND stderr to the single PTY
// slave, so the stderr TTY check passes too.
//
// When stdout is not a TTY (CI, piped, --print-cmd consumers) we keep plain
// inherited fds, so non-interactive behavior is unchanged. The run_analysis.sh
// argument vector is byte-identical on both paths (parity with
// run_interactive.sh --print-cmd is a maintained contract).
//
// POSIX-only (creack/pty + SIGWINCH); the project targets macOS + Linux.
func RunPipeline(projectDir string, c Config) error {
	fmt.Println(mutedStyle.Render("Launching… (full pipeline output follows)"))
	fmt.Println(mutedStyle.Render(c.CommandPreview()))
	fmt.Println()

	cmd := exec.Command("bash", append([]string{"run_analysis.sh"}, c.ToArgs()...)...)
	cmd.Dir = projectDir

	// Plain inherited fds unless stdout is an interactive terminal.
	if !isatty.IsTerminal(os.Stdout.Fd()) {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		cmd.Stdin = os.Stdin
		return cmd.Run()
	}

	// Interactive terminal: run under a PTY so cli renders live progress.
	ptmx, err := pty.Start(cmd)
	if err != nil {
		// PTY allocation failed — degrade gracefully to plain fds.
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		cmd.Stdin = os.Stdin
		return cmd.Run()
	}
	defer func() { _ = ptmx.Close() }()

	// Keep the PTY sized to the real terminal (initial + on window resize).
	winch := make(chan os.Signal, 1)
	signal.Notify(winch, syscall.SIGWINCH)
	go func() {
		for range winch {
			_ = pty.InheritSize(os.Stdout, ptmx)
		}
	}()
	winch <- syscall.SIGWINCH // trigger the initial sizing
	defer func() { signal.Stop(winch); close(winch) }()

	// Best-effort forward stdin (the pipeline reads none mid-run, but a
	// launcher confirm prompt might); mirror child output to our stdout.
	go func() { _, _ = io.Copy(ptmx, os.Stdin) }()
	_, _ = io.Copy(os.Stdout, ptmx)

	return cmd.Wait()
}
