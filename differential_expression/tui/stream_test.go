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

func TestFlushAllDrainsEverything(t *testing.T) {
	m := runModel{
		queue:   []string{"a", "b", "c"},
		cur:     []rune("in flight"),
		curKind: kindInfo,
		pos:     3,
	}
	cmd := m.flushAll()
	if cmd == nil {
		t.Fatal("flushAll returned no command")
	}
	if len(m.queue) != 0 || len(m.cur) != 0 || m.pos != 0 {
		t.Errorf("state not drained: queue=%d cur=%d pos=%d", len(m.queue), len(m.cur), m.pos)
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
