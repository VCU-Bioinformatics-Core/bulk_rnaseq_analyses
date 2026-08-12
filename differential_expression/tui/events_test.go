package main

import (
	"context"
	"os"
	"testing"
	"time"
)

func TestParseEvent(t *testing.T) {
	cases := []struct {
		name string
		line string
		want Event
	}{
		{"start", `{"t":"start","total":30,"comparisons":3,"runid":"r1"}`,
			Event{T: "start", Total: 30, Comparisons: 3, RunID: "r1"}},
		{"tick", `{"t":"tick","current":7,"total":30,"i":1,"n":3,"comparison":"A_vs_B","step":"volcano"}`,
			Event{T: "tick", Current: 7, Total: 30, I: 1, N: 3, Comparison: "A_vs_B", Step: "volcano"}},
		{"phase", `{"t":"phase","name":"PCA"}`, Event{T: "phase", Name: "PCA"}},
		{"done", `{"t":"done","ok":true}`, Event{T: "done", OK: true}},
		{"trailing whitespace tolerated", "  {\"t\":\"done\",\"ok\":true}  \n", Event{T: "done", OK: true}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := ParseEvent(tc.line)
			if err != nil {
				t.Fatalf("ParseEvent(%q) error: %v", tc.line, err)
			}
			if got != tc.want {
				t.Errorf("got %+v, want %+v", got, tc.want)
			}
		})
	}
}

func TestParseEventRejectsMalformed(t *testing.T) {
	if _, err := ParseEvent(`{"t":"tick",`); err == nil {
		t.Error("expected an error for truncated JSON")
	}
}

func TestTailEventsFollowsAppendsAndSkipsGarbage(t *testing.T) {
	f, err := os.CreateTemp(t.TempDir(), "events_*.ndjson")
	if err != nil {
		t.Fatal(err)
	}
	path := f.Name()

	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	ch := make(chan Event, 16)
	go TailEvents(ctx, path, ch)

	// Written after the tailer starts: it must follow appends, not just read once.
	writes := []string{
		`{"t":"start","total":20,"comparisons":2}` + "\n",
		"not json at all\n", // must be skipped, not fatal
		`{"t":"tick","current":5,"total":20,"step":"DESeq2"}` + "\n",
		`{"t":"done","ok":true}` + "\n",
	}
	go func() {
		for _, w := range writes {
			time.Sleep(20 * time.Millisecond)
			_, _ = f.WriteString(w)
		}
		_ = f.Close()
	}()

	var got []Event
	deadline := time.After(5 * time.Second)
	for len(got) < 3 {
		select {
		case e := <-ch:
			got = append(got, e)
		case <-deadline:
			t.Fatalf("timed out; got %d events: %+v", len(got), got)
		}
	}

	if got[0].T != "start" || got[0].Total != 20 {
		t.Errorf("first event = %+v", got[0])
	}
	if got[1].T != "tick" || got[1].Current != 5 || got[1].Step != "DESeq2" {
		t.Errorf("second event (garbage line should have been skipped) = %+v", got[1])
	}
	if got[2].T != "done" || !got[2].OK {
		t.Errorf("third event = %+v", got[2])
	}
}

func TestTailEventsWaitsForAMissingFile(t *testing.T) {
	// The tailer may start before the writer creates the file.
	path := t.TempDir() + "/not_yet.ndjson"
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	ch := make(chan Event, 4)
	go TailEvents(ctx, path, ch)

	time.Sleep(120 * time.Millisecond)
	if err := os.WriteFile(path, []byte(`{"t":"phase","name":"PCA"}`+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	select {
	case e := <-ch:
		if e.T != "phase" || e.Name != "PCA" {
			t.Errorf("got %+v", e)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("tailer never picked up the late-created file")
	}
}

func TestTailEventsStopsOnContextCancel(t *testing.T) {
	path := t.TempDir() + "/ev.ndjson"
	if err := os.WriteFile(path, []byte(""), 0o644); err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	ch := make(chan Event)
	go TailEvents(ctx, path, ch)
	cancel()

	select {
	case _, open := <-ch:
		if open {
			t.Error("expected the channel to be closed after cancel")
		}
	case <-time.After(3 * time.Second):
		t.Fatal("tailer did not stop on context cancel")
	}
}
