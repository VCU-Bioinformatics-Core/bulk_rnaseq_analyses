package main

import (
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/charmbracelet/lipgloss"
)

// Output styling + typewriter streaming.
//
// Two jobs, and they are related:
//
//  1. Colour. Pipeline lines are classified by their `cli` glyph and given a
//     deliberate, consistent style, so severity is readable at a glance
//     instead of being a wall of one colour.
//
//  2. Streaming. Lines are typed out a few characters at a time rather than
//     dumped whole. This REQUIRES stripping the incoming ANSI first: emitting
//     an escape sequence one byte at a time would split it and spray garbage
//     across the terminal. So we strip R's colours and re-apply our own — which
//     is exactly what makes (1) possible.

var ansiRE = regexp.MustCompile(`\x1b\[[0-9;?]*[ -/]*[@-~]`)

// stripANSI removes CSI escape sequences so a line can be safely typed one
// character at a time and re-styled by us.
func stripANSI(s string) string { return ansiRE.ReplaceAllString(s, "") }

type lineKind int

const (
	kindPlain lineKind = iota
	kindInfo
	kindSuccess
	kindWarn
	kindError
	kindHeader
)

var (
	styPlain   = lipgloss.NewStyle()
	styInfo    = lipgloss.NewStyle().Foreground(lipgloss.Color("111"))
	stySuccess = lipgloss.NewStyle().Foreground(lipgloss.Color("78"))
	styWarn    = lipgloss.NewStyle().Foreground(lipgloss.Color("214"))
	styError   = lipgloss.NewStyle().Bold(true).Foreground(lipgloss.Color("203"))
	styHeader  = lipgloss.NewStyle().Bold(true).Foreground(lipgloss.Color("147"))
)

// classifyLine maps a pipeline output line to a severity. It matches the
// glyphs cli emits in both its unicode and ASCII fallback modes.
func classifyLine(s string) lineKind {
	t := strings.TrimSpace(stripANSI(s))
	if t == "" {
		return kindPlain
	}
	switch {
	case strings.HasPrefix(t, "──"), strings.HasPrefix(t, "--- "), strings.HasPrefix(t, "== "):
		return kindHeader
	case strings.HasPrefix(t, "✔"), strings.HasPrefix(t, "v "):
		return kindSuccess
	case strings.HasPrefix(t, "ℹ"), strings.HasPrefix(t, "i "):
		return kindInfo
	case strings.HasPrefix(t, "!"), strings.HasPrefix(t, "⚠"):
		return kindWarn
	case strings.HasPrefix(t, "✖"), strings.HasPrefix(t, "✗"), strings.HasPrefix(t, "x "),
		strings.HasPrefix(t, "Error"), strings.HasPrefix(t, "Fatal"):
		return kindError
	}
	return kindPlain
}

func styleFor(k lineKind) lipgloss.Style {
	switch k {
	case kindInfo:
		return styInfo
	case kindSuccess:
		return stySuccess
	case kindWarn:
		return styWarn
	case kindError:
		return styError
	case kindHeader:
		return styHeader
	}
	return styPlain
}

// typeDelayMS is the per-character delay for streamed output. Configurable via
// BISR_TUI_TYPE_DELAY_MS; 0 disables typing and prints lines immediately.
func typeDelayMS() int {
	if v, ok := os.LookupEnv("BISR_TUI_TYPE_DELAY_MS"); ok {
		if n, err := strconv.Atoi(strings.TrimSpace(v)); err == nil && n >= 0 {
			return n
		}
	}
	return 12
}

// typeChunk decides how many characters to reveal per tick.
//
// A fixed 1-char-per-tick reads beautifully for a trickle of status lines but
// falls catastrophically behind a burst (Quarto alone dumps hundreds of lines),
// which would leave the UI narrating output long after the run finished. So the
// reveal rate scales with the backlog, and past a threshold we stop typing
// altogether and flush — the animation must never delay the pipeline's story.
func typeChunk(backlog int) int {
	switch {
	case backlog >= flushBacklog:
		return 0 // caller flushes the whole queue instantly
	case backlog > 8:
		return 8
	case backlog > 3:
		return 3
	}
	return 1
}

// flushBacklog is the queue length past which typing is abandoned for this
// burst and everything pending is printed at once.
const flushBacklog = 40
