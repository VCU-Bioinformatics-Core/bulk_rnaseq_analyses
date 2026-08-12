package main

import (
	"fmt"
	"strings"
	"time"

	"github.com/charmbracelet/bubbles/progress"
	"github.com/charmbracelet/bubbles/spinner"
	tea "github.com/charmbracelet/bubbletea"
	"github.com/charmbracelet/lipgloss"
)

// Live run UI.
//
// The pipeline's structured events (events.go) drive a bubbletea model, so the
// bar is redrawn from STATE rather than scraped from R's pretty output — the
// same "state -> re-render" idea Claude Code gets from React/Ink. The child's
// own log lines are printed ABOVE the persistent UI with tea.Println, so the
// progress block stays pinned to the bottom while output scrolls past it.
//
// Messages arriving from the two goroutines started by runner.go:
type eventMsg Event      // one decoded NDJSON progress event
type outputMsg string    // one line of child stdout/stderr
type streamDoneMsg struct{ which string } // a stream closed
type finishedMsg struct{ err error }      // the child process exited
type typeTickMsg struct{}                 // reveal the next characters

var (
	stepStyle  = lipgloss.NewStyle().Foreground(lipgloss.Color("246"))
	nameStyle  = lipgloss.NewStyle().Bold(true).Foreground(lipgloss.Color("39"))
	phaseStyle = lipgloss.NewStyle().Foreground(lipgloss.Color("213"))
	failStyle  = lipgloss.NewStyle().Bold(true).Foreground(lipgloss.Color("196"))
)

type runModel struct {
	prog progress.Model
	spin spinner.Model

	events <-chan Event
	output <-chan string
	done   <-chan error

	total, current int
	comparison     string
	step           string
	phase          string
	i, n           int

	// Typewriter streaming: `cur` is the line being revealed, `queue` the
	// lines waiting behind it. Completed lines are committed to scrollback
	// with tea.Println; only the in-flight one lives in the View.
	queue   []string
	cur     []rune
	curKind lineKind
	pos     int
	delay   time.Duration

	finished bool
	err      error
	// Once both the child has exited and the event stream has closed there is
	// nothing left to render, so the program can quit.
	childExited  bool
	eventsClosed bool
}

func newRunModel(events <-chan Event, output <-chan string, done <-chan error) runModel {
	p := progress.New(progress.WithDefaultGradient(), progress.WithWidth(38))
	s := spinner.New()
	s.Spinner = spinner.MiniDot
	s.Style = lipgloss.NewStyle().Foreground(lipgloss.Color("39"))
	return runModel{
		prog: p, spin: s,
		events: events, output: output, done: done,
		phase: "starting", step: "…",
		delay: time.Duration(typeDelayMS()) * time.Millisecond,
	}
}

// printLine commits one finished line to scrollback, styled by severity.
func printLine(text string, kind lineKind) tea.Cmd {
	return tea.Println(styleFor(kind).Render(text))
}

// startNext begins revealing the next queued line, or reports that the
// typewriter is now idle.
func (m *runModel) startNext() tea.Cmd {
	if len(m.queue) == 0 {
		m.cur, m.pos = nil, 0
		return nil
	}
	next := m.queue[0]
	m.queue = m.queue[1:]
	m.cur = []rune(next)
	m.curKind = classifyLine(next)
	m.pos = 0
	return tea.Tick(m.delay, func(time.Time) tea.Msg { return typeTickMsg{} })
}

// flushAll abandons typing for this burst and prints everything pending at
// once, so a flood of output can never make the UI lag behind the pipeline.
func (m *runModel) flushAll() tea.Cmd {
	cmds := make([]tea.Cmd, 0, len(m.queue)+1)
	if len(m.cur) > 0 {
		cmds = append(cmds, printLine(string(m.cur), m.curKind))
	}
	for _, l := range m.queue {
		cmds = append(cmds, printLine(l, classifyLine(l)))
	}
	m.queue, m.cur, m.pos = nil, nil, 0
	return tea.Batch(cmds...)
}

func (m runModel) Init() tea.Cmd {
	return tea.Batch(m.spin.Tick, waitEvent(m.events), waitOutput(m.output), waitDone(m.done))
}

func waitEvent(ch <-chan Event) tea.Cmd {
	return func() tea.Msg {
		e, ok := <-ch
		if !ok {
			return streamDoneMsg{which: "events"}
		}
		return eventMsg(e)
	}
}

func waitOutput(ch <-chan string) tea.Cmd {
	return func() tea.Msg {
		s, ok := <-ch
		if !ok {
			return streamDoneMsg{which: "output"}
		}
		return outputMsg(s)
	}
}

func waitDone(ch <-chan error) tea.Cmd {
	return func() tea.Msg {
		err, ok := <-ch
		if !ok {
			return finishedMsg{}
		}
		return finishedMsg{err: err}
	}
}

func (m runModel) Update(msg tea.Msg) (tea.Model, tea.Cmd) {
	switch msg := msg.(type) {

	case eventMsg:
		switch msg.T {
		case "start":
			m.total = msg.Total
			m.n = msg.Comparisons
			m.phase = "differential expression"
		case "tick":
			m.current, m.step = msg.Current, msg.Step
			m.comparison, m.i = msg.Comparison, msg.I
			if msg.Total > 0 {
				m.total = msg.Total
			}
			if msg.N > 0 {
				m.n = msg.N
			}
		case "phase":
			m.phase, m.step = msg.Name, "…"
		case "done":
			m.finished = true
			if m.total > 0 {
				m.current = m.total
			}
		case "error":
			m.err = fmt.Errorf("%s", msg.Msg)
		}
		return m, waitEvent(m.events)

	case outputMsg:
		// Strip R's ANSI so the line can be typed a character at a time without
		// splitting an escape sequence; we re-colour it by severity ourselves.
		line := stripANSI(strings.TrimRight(string(msg), "\r\n"))
		next := waitOutput(m.output)

		if m.delay <= 0 { // typing disabled: straight to scrollback
			return m, tea.Batch(printLine(line, classifyLine(line)), next)
		}
		m.queue = append(m.queue, line)
		if len(m.queue) >= flushBacklog {
			return m, tea.Batch(m.flushAll(), next)
		}
		if len(m.cur) == 0 { // idle -> start revealing immediately
			return m, tea.Batch(m.startNext(), next)
		}
		return m, next

	case typeTickMsg:
		if len(m.cur) == 0 {
			return m, m.startNext()
		}
		chunk := typeChunk(len(m.queue))
		if chunk == 0 { // backlog blew past the threshold mid-line
			return m, m.flushAll()
		}
		m.pos += chunk
		if m.pos < len(m.cur) {
			return m, tea.Tick(m.delay, func(time.Time) tea.Msg { return typeTickMsg{} })
		}
		// Line finished: commit it to scrollback and move on.
		done := printLine(string(m.cur), m.curKind)
		return m, tea.Batch(done, m.startNext())

	case streamDoneMsg:
		if msg.which == "events" {
			m.eventsClosed = true
		}
		if m.childExited && m.eventsClosed {
			return m, tea.Quit
		}
		return m, nil

	case finishedMsg:
		m.childExited = true
		if msg.err != nil {
			m.err = msg.err
		}
		// Don't drop anything still being typed when the child exits.
		return m, tea.Sequence(m.flushAll(), tea.Quit)

	case spinner.TickMsg:
		var cmd tea.Cmd
		m.spin, cmd = m.spin.Update(msg)
		return m, cmd

	case tea.KeyMsg:
		// Ctrl-C is forwarded to the child by the terminal (it shares our
		// process group), so just let the UI wind down.
		if msg.Type == tea.KeyCtrlC {
			return m, tea.Quit
		}
	}
	return m, nil
}

func (m runModel) View() string {
	// Nothing lingering once the run is over — the final summary is printed by
	// runner.go after the program exits.
	if m.childExited {
		return ""
	}

	pct := 0.0
	if m.total > 0 {
		pct = float64(m.current) / float64(m.total)
		if pct > 1 {
			pct = 1
		}
	}

	head := phaseStyle.Render(m.phase)
	if m.comparison != "" {
		head = fmt.Sprintf("%s %s", nameStyle.Render(m.comparison),
			stepStyle.Render(fmt.Sprintf("(%d/%d)", m.i, m.n)))
	}

	var b strings.Builder
	// The line currently being revealed lives here until it is complete, then
	// it is committed to scrollback with tea.Println.
	if m.pos > 0 && m.pos < len(m.cur) {
		b.WriteString(styleFor(m.curKind).Render(string(m.cur[:m.pos])) + "\n")
	}
	b.WriteString("\n")
	b.WriteString(fmt.Sprintf("%s %s  %s\n", m.spin.View(), head, stepStyle.Render("· "+m.step)))
	b.WriteString("  " + m.prog.ViewAs(pct))
	if m.total > 0 {
		b.WriteString(stepStyle.Render(fmt.Sprintf("  %d/%d", m.current, m.total)))
	}
	if m.err != nil {
		b.WriteString("\n  " + failStyle.Render("✗ "+m.err.Error()))
	}
	b.WriteString("\n")
	return b.String()
}
