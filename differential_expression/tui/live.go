package main

import (
	"bufio"
	"context"
	"errors"
	"os"
	"os/exec"
	"os/signal"
	"sync"
	"syscall"

	tea "github.com/charmbracelet/bubbletea"
	"github.com/creack/pty"
)

// ErrLiveUIUnavailable signals that the event-driven UI could not be started
// and the caller should fall back to plainly relaying the child's output. The
// analysis must never fail because the UI did.
var ErrLiveUIUnavailable = errors.New("live UI unavailable")

// runLiveUI starts cmd under a pseudo-terminal and renders a bubbletea UI fed
// by two streams: the pipeline's NDJSON progress events (tailed from
// eventsPath) and the child's own output, which is printed above the pinned
// progress block. Returns the child's exit error.
func runLiveUI(cmd *exec.Cmd, eventsPath string) error {
	ptmx, err := pty.Start(cmd)
	if err != nil {
		return ErrLiveUIUnavailable
	}
	defer func() { _ = ptmx.Close() }()

	// Keep the child's terminal the same size as ours.
	winch := make(chan os.Signal, 1)
	signal.Notify(winch, syscall.SIGWINCH)
	go func() {
		for range winch {
			_ = pty.InheritSize(os.Stdout, ptmx)
		}
	}()
	winch <- syscall.SIGWINCH
	defer func() { signal.Stop(winch); close(winch) }()

	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	// Progress events.
	events := make(chan Event, 64)
	go TailEvents(ctx, eventsPath, events)

	// Child output, line by line.
	output := make(chan string, 256)
	go func() {
		defer close(output)
		sc := bufio.NewScanner(ptmx)
		sc.Buffer(make([]byte, 0, 64*1024), 1024*1024)
		for sc.Scan() {
			select {
			case <-ctx.Done():
				return
			case output <- sc.Text():
			}
		}
	}()

	// Child exit.
	var (
		wg     sync.WaitGroup
		runErr error
	)
	done := make(chan error, 1)
	wg.Add(1)
	go func() {
		defer wg.Done()
		e := cmd.Wait()
		runErr = e
		done <- e
		close(done)
	}()

	p := tea.NewProgram(newRunModel(events, output, done))
	if _, perr := p.Run(); perr != nil {
		// The UI died; still let the child finish so we report its real status.
		wg.Wait()
		if runErr != nil {
			return runErr
		}
		return perr
	}

	cancel()
	wg.Wait() // establishes visibility of runErr
	return runErr
}
