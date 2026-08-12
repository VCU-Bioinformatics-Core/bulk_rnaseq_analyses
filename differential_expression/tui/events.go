package main

import (
	"bufio"
	"context"
	"encoding/json"
	"io"
	"os"
	"strings"
	"time"
)

// Event is one line of the pipeline's NDJSON progress stream. The R side
// (bisrDE/R/events.R) appends these when BISR_EVENTS_FILE is set; fields not
// relevant to a given event type are simply absent and stay zero-valued.
//
//	{"t":"start","total":30,"comparisons":3,"runid":"..."}
//	{"t":"phase","name":"PCA"}
//	{"t":"tick","current":7,"total":30,"i":1,"n":3,"comparison":"A_vs_B","step":"volcano"}
//	{"t":"done","ok":true}
type Event struct {
	T           string `json:"t"`
	Total       int    `json:"total"`
	Current     int    `json:"current"`
	Comparisons int    `json:"comparisons"`
	I           int    `json:"i"`
	N           int    `json:"n"`
	Comparison  string `json:"comparison"`
	Step        string `json:"step"`
	Name        string `json:"name"`
	RunID       string `json:"runid"`
	OK          bool   `json:"ok"`
	Msg         string `json:"msg"`
}

// ParseEvent decodes a single NDJSON line. Unparseable lines are reported as
// an error so the caller can skip them rather than abort the UI — a malformed
// telemetry line must never kill a running analysis.
func ParseEvent(line string) (Event, error) {
	var e Event
	err := json.Unmarshal([]byte(strings.TrimSpace(line)), &e)
	return e, err
}

// TailEvents follows an NDJSON file, emitting each decoded Event on out until
// ctx is cancelled. It tolerates the file not existing yet, and buffers
// partial trailing lines so a half-flushed write is never parsed. The channel
// is closed on return.
func TailEvents(ctx context.Context, path string, out chan<- Event) {
	defer close(out)

	var f *os.File
	// The writer may not have created the file yet when we start.
	for f == nil {
		if fh, err := os.Open(path); err == nil {
			f = fh
		} else {
			select {
			case <-ctx.Done():
				return
			case <-time.After(50 * time.Millisecond):
			}
		}
	}
	defer f.Close()

	reader := bufio.NewReader(f)
	var partial strings.Builder

	for {
		line, err := reader.ReadString('\n')
		if err == io.EOF {
			if line != "" {
				// Incomplete final line: hold it until the newline arrives.
				partial.WriteString(line)
			}
			select {
			case <-ctx.Done():
				return
			case <-time.After(60 * time.Millisecond):
				continue
			}
		}
		if err != nil {
			return
		}
		if partial.Len() > 0 {
			line = partial.String() + line
			partial.Reset()
		}
		if strings.TrimSpace(line) == "" {
			continue
		}
		ev, perr := ParseEvent(line)
		if perr != nil {
			continue // skip malformed telemetry, keep the UI alive
		}
		select {
		case <-ctx.Done():
			return
		case out <- ev:
		}
	}
}
