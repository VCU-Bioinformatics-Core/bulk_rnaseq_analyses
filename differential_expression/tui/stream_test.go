package main

import (
	"os"
	"strings"
	"testing"
)

func TestStripANSI(t *testing.T) {
	cases := map[string]string{
		"\x1b[32mv Loaded 10 genes\x1b[39m": "v Loaded 10 genes",
		"\x1b[1m\x1b[36m-- PCA --\x1b[0m":   "-- PCA --",
		"plain text":                        "plain text",
		"":                                  "",
	}
	for in, want := range cases {
		if got := stripANSI(in); got != want {
			t.Errorf("stripANSI(%q) = %q, want %q", in, got, want)
		}
	}
}

func TestStripANSILeavesNoEscapeBytes(t *testing.T) {
	// Critical for the typewriter: any surviving escape byte could be split
	// mid-sequence when revealing characters one at a time.
	in := "\x1b[33m! warning\x1b[39m \x1b[1mbold\x1b[0m"
	if got := stripANSI(in); strings.ContainsRune(got, '\x1b') {
		t.Errorf("escape byte survived: %q", got)
	}
}

func TestClassifyLine(t *testing.T) {
	cases := []struct {
		line string
		want lineKind
	}{
		{"── Loading data ──────────", kindHeader},
		{"v Loaded 10 genes x 60 samples", kindSuccess},
		{"✔ RDS saved", kindSuccess},
		{"i Annotating results...", kindInfo},
		{"ℹ Session log: /tmp/x.log", kindInfo},
		{"! No DisplayName column", kindWarn},
		{"✖ run_analysis error", kindError},
		{"Error in foo(): bad", kindError},
		{"just some output", kindPlain},
		{"", kindPlain},
		// Styling must survive R's colours.
		{"\x1b[32mv Loaded\x1b[39m", kindSuccess},
	}
	for _, tc := range cases {
		if got := classifyLine(tc.line); got != tc.want {
			t.Errorf("classifyLine(%q) = %v, want %v", tc.line, got, tc.want)
		}
	}
}

func TestTypeChunkScalesWithBacklog(t *testing.T) {
	// Small backlog: one character at a time, for the organic feel.
	if got := typeChunk(0); got != 1 {
		t.Errorf("typeChunk(0) = %d, want 1", got)
	}
	if got := typeChunk(2); got != 1 {
		t.Errorf("typeChunk(2) = %d, want 1", got)
	}
	// Growing backlog: accelerate so the UI keeps up.
	if got := typeChunk(5); got <= 1 {
		t.Errorf("typeChunk(5) = %d, want > 1", got)
	}
	if got := typeChunk(20); got < typeChunk(5) {
		t.Errorf("reveal rate must not decrease as backlog grows")
	}
	// Flood: stop typing entirely and flush.
	if got := typeChunk(flushBacklog); got != 0 {
		t.Errorf("typeChunk(%d) = %d, want 0 (flush)", flushBacklog, got)
	}
}

func TestTypeDelayMSConfigurable(t *testing.T) {
	t.Setenv("BISR_TUI_TYPE_DELAY_MS", "25")
	if got := typeDelayMS(); got != 25 {
		t.Errorf("typeDelayMS() = %d, want 25", got)
	}
	t.Setenv("BISR_TUI_TYPE_DELAY_MS", "0") // disables typing
	if got := typeDelayMS(); got != 0 {
		t.Errorf("typeDelayMS() = %d, want 0", got)
	}
	t.Setenv("BISR_TUI_TYPE_DELAY_MS", "garbage") // falls back to default
	if got := typeDelayMS(); got != 12 {
		t.Errorf("typeDelayMS() = %d, want default 12", got)
	}
	os.Unsetenv("BISR_TUI_TYPE_DELAY_MS")
	if got := typeDelayMS(); got != 12 {
		t.Errorf("typeDelayMS() unset = %d, want default 12", got)
	}
}

func TestFlushAllDrainsIntoTheRing(t *testing.T) {
	m := runModel{
		queue:    []string{"a", "b", "c"},
		cur:      []rune("in flight"),
		curKind:  kindInfo,
		pos:      3,
		maxLines: 8,
	}
	m.flushAll()

	if len(m.queue) != 0 || len(m.cur) != 0 || m.pos != 0 {
		t.Errorf("state not drained: queue=%d cur=%d pos=%d", len(m.queue), len(m.cur), m.pos)
	}
	// The in-flight line plus all three queued lines must survive.
	if len(m.ring) != 4 {
		t.Fatalf("ring has %d lines, want 4: %+v", len(m.ring), m.ring)
	}
	if m.ring[0].text != "in flight" || m.ring[3].text != "c" {
		t.Errorf("unexpected ring contents: %+v", m.ring)
	}
}

func TestRingBufferIsBounded(t *testing.T) {
	m := runModel{maxLines: 5}
	for i := 0; i < 20; i++ {
		m.push(string(rune('a'+i)), kindPlain)
	}
	if len(m.ring) != 5 {
		t.Fatalf("ring grew to %d, want a hard cap of 5", len(m.ring))
	}
	// It must keep the NEWEST lines, dropping the oldest off the top.
	if m.ring[0].text != "p" || m.ring[4].text != "t" {
		t.Errorf("ring kept the wrong window: %+v", m.ring)
	}
}

func TestBoxHasStableHeight(t *testing.T) {
	m := runModel{maxLines: 8}
	if got := len(m.boxLines(60)); got != 8 {
		t.Errorf("empty box rendered %d rows, want 8 (padded)", got)
	}
	for i := 0; i < 3; i++ {
		m.push("line", kindInfo)
	}
	if got := len(m.boxLines(60)); got != 8 {
		t.Errorf("partly filled box rendered %d rows, want 8", got)
	}
	for i := 0; i < 50; i++ {
		m.push("line", kindInfo)
	}
	if got := len(m.boxLines(60)); got != 8 {
		t.Errorf("overfull box rendered %d rows, want 8", got)
	}
}

func TestTruncate(t *testing.T) {
	if got := truncate("hello", 10); got != "hello" {
		t.Errorf("short line altered: %q", got)
	}
	got := truncate("abcdefghij", 5)
	if []rune(got)[len([]rune(got))-1] != '…' {
		t.Errorf("truncated line should end in an ellipsis, got %q", got)
	}
	if len([]rune(got)) != 5 {
		t.Errorf("truncate(_, 5) produced %d runes: %q", len([]rune(got)), got)
	}
	if truncate("anything", 1) != "" {
		t.Error("degenerate width should produce an empty string, not a panic")
	}
	// Multi-byte content must be cut by rune, never mid-character.
	if got := truncate("αβγδεζηθ", 4); len([]rune(got)) != 4 {
		t.Errorf("unicode truncate produced %d runes: %q", len([]rune(got)), got)
	}
}

func TestTUILinesClamped(t *testing.T) {
	t.Setenv("BISR_TUI_LINES", "10")
	if got := tuiLines(); got != 10 {
		t.Errorf("tuiLines() = %d, want 10", got)
	}
	t.Setenv("BISR_TUI_LINES", "1") // too small to be useful
	if got := tuiLines(); got != 3 {
		t.Errorf("tuiLines() = %d, want clamp to 3", got)
	}
	t.Setenv("BISR_TUI_LINES", "9999") // would fill the screen
	if got := tuiLines(); got != 40 {
		t.Errorf("tuiLines() = %d, want clamp to 40", got)
	}
	t.Setenv("BISR_TUI_LINES", "nonsense")
	if got := tuiLines(); got != 8 {
		t.Errorf("tuiLines() = %d, want default 8", got)
	}
}

func TestStartNextAdvancesTheQueue(t *testing.T) {
	m := runModel{queue: []string{"v first", "i second"}}
	if cmd := m.startNext(); cmd == nil {
		t.Fatal("expected a tick command")
	}
	if string(m.cur) != "v first" {
		t.Errorf("cur = %q, want %q", string(m.cur), "v first")
	}
	if m.curKind != kindSuccess {
		t.Errorf("curKind = %v, want kindSuccess", m.curKind)
	}
	if len(m.queue) != 1 {
		t.Errorf("queue not advanced: %v", m.queue)
	}
	// Draining to empty must report idle (nil command), not spin forever.
	m.queue = nil
	m.startNext()
	if cmd := m.startNext(); cmd != nil {
		t.Error("expected nil command when the queue is empty")
	}
}
