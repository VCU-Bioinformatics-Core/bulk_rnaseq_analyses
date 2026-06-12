# bisrde-tui — Go TUI launcher (charmbracelet Phase B)

A [bubbletea](https://github.com/charmbracelet/bubbletea)/[huh](https://github.com/charmbracelet/huh)/[lipgloss](https://github.com/charmbracelet/lipgloss)
front-end for the bisrDE bulk RNA-seq DE pipeline. It collects parameters
through a styled form — including interactive **group / sample / contrast**
selection derived from the samplesheet — then hands off to `run_analysis.sh`
with the **exact same flag vector** that `run_interactive.sh` assembles.

It is the "Go TUI" option offered by `run_interactive.sh`; you can also run it
directly.

## Build

Needs Go 1.23+ (`brew install go`):

```bash
cd tui
go build -o bisrde-tui .
```

Cross-compile a Linux x86_64 binary for HPC (from a mac):

```bash
GOOS=linux GOARCH=amd64 go build -o bisrde-tui-linux-amd64 .
```

The binary is git-ignored (built per-platform). Ship prebuilt binaries as
GitHub release assets so HPC users don't need Go; the container
(`dge_analysis.def`) can `COPY` the Linux binary in.

## Run

```bash
# via the chooser:
bash ../run_interactive.sh        # pick "Go TUI"

# or directly (from the differential_expression/ project dir):
./tui/bisrde-tui

# assemble + print the run_analysis.sh command without running (for tests):
./tui/bisrde-tui --print-cmd
```

Set `BISR_PROJECT_DIR` to the `differential_expression/` dir if you launch the
binary from elsewhere (the chooser does this for you). `BISR_*` env vars
pre-answer prompts (same names as `run_interactive.sh`).

## Design notes

- **Form, not a fake terminal.** The pipeline emits rich `cli` progress
  (banners, progress bars, colour). The TUI's value is the styled huh form +
  lipgloss panels around the run; it streams the pipeline's own output live
  and unmodified rather than re-rendering it. A live multi-pane log viewport
  (bubbles `viewport`) is a possible future enhancement.
- **Flag parity is a contract.** `Config.ToArgs()` (config.go) must produce
  the same `run_analysis.sh` argument vector as `run_interactive.sh
  --print-cmd`. There is a parity check in the project tracker; keep the two
  in sync when adding flags.
- **Metadata columns** recognised by `samplesheet.go` (`SampleID`, `GroupID`,
  `Exclude`, `DisplayName`) mirror bisrDE's `.samplesheet_meta_cols` in
  `bisrDE/R/io.R`.

## Files

| File             | Role |
|------------------|------|
| `main.go`        | entrypoint: flags, TTY detection, form → confirm → run |
| `form.go`        | huh forms + `BISR_*` env overrides → `Config` |
| `config.go`      | `Config` struct, `ToArgs()` flag assembly, summary |
| `samplesheet.go` | CSV parse → groups / samples / contrasts |
| `runner.go`      | `exec` `run_analysis.sh`, native streaming output |
| `theme.go`       | lipgloss palette (matches `run_interactive.sh`) |
