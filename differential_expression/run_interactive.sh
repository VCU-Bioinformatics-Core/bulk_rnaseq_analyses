#!/usr/bin/env bash
#
# run_interactive.sh — beautiful interactive front-door for the bisrDE
# bulk RNA-seq DE pipeline (v1.5.2).
#
# Uses charmbracelet tools when available and degrades gracefully to plain
# shell prompts otherwise:
#   - gum   : styled banner/borders (lipgloss), forms (huh), spinner +
#             multi-select (bubbles). brew install gum
#   - glow  : renders the run summary as terminal markdown. brew install glow
#   - freeze: (docs) screenshots the CLI. brew install freeze
#
# It collects every parameter — including the v1.5.0 sample/group/contrast
# selectors — then hands off to ./run_analysis.sh with the assembled flags.
#
# Testing / scripting hooks (no TTY, no gum needed):
#   - Set BISR_* env vars to pre-answer prompts (skips the interactive ask).
#   - Pass --print-cmd to print the assembled run_analysis.sh command and
#     exit instead of executing it.
#
# Usage:
#   bash run_interactive.sh                 # full interactive experience
#   bash run_interactive.sh --print-cmd     # assemble + print (for testing)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Theme (256-color codes; override via env).
# ---------------------------------------------------------------------------
BISR_PRIMARY="${BISR_PRIMARY:-212}"   # pink
BISR_ACCENT="${BISR_ACCENT:-86}"      # cyan
BISR_MUTED="${BISR_MUTED:-244}"       # grey

PRINT_CMD=0
for arg in "$@"; do
  case "$arg" in
    --print-cmd) PRINT_CMD=1 ;;
    -h|--help)
      cat <<'EOF'
run_interactive.sh — interactive launcher for the bisrDE pipeline.

Run it with no arguments for the full guided experience. It collects the
counts file, samplesheet, annotation, run ID, output dir, optional BRS
ticket and ID type, then offers interactive sample / group / contrast
exclusion before launching the analysis.

Install the (optional) charmbracelet tools for the styled experience:
    brew install gum glow freeze

Without them, the launcher falls back to plain prompts.

Flags:
    --print-cmd   Assemble the run_analysis.sh command and print it (no run).
    -h, --help    This help.
EOF
      exit 0 ;;
  esac
done

have() { command -v "$1" >/dev/null 2>&1; }

# Interactive only when we have a real terminal and aren't just printing the
# command. In non-interactive mode the ask_* helpers return presets/defaults
# instead of blocking on read/gum (keeps the script scriptable + testable).
INTERACTIVE=1
{ [ "$PRINT_CMD" = "1" ] || [ ! -t 0 ]; } && INTERACTIVE=0

# ---------------------------------------------------------------------------
# Prompt primitives (gum when present, plain shell otherwise).
# Each respects a pre-set value passed as $3 (from an env override) so the
# whole flow is scriptable/testable without a TTY.
# ---------------------------------------------------------------------------
banner() {
  if have gum; then
    gum style --border double --margin "1 0" --padding "1 4" \
      --border-foreground "$BISR_PRIMARY" --foreground "$BISR_PRIMARY" \
      "Bulk RNA-Seq Differential Expression" "bisrDE · v1.5.2 · VCU Massey BISR"
  else
    echo "======================================================="
    echo "  Bulk RNA-Seq DE Pipeline — bisrDE v1.5.2 (VCU Massey)"
    echo "======================================================="
  fi
}

note() {  # styled muted line
  if have gum; then gum style --foreground "$BISR_MUTED" "$1"; else echo "$1"; fi
}

ask_input() {  # prompt, default, preset -> echoes value
  local prompt="$1" default="${2:-}" preset="${3:-}"
  if [ -n "$preset" ]; then echo "$preset"; return 0; fi
  if [ "$INTERACTIVE" = "0" ]; then echo "$default"; return 0; fi
  if have gum; then
    gum input --prompt "$prompt › " --placeholder "$default" --value "$default"
  else
    local v; read -r -p "$prompt [$default]: " v; echo "${v:-$default}"
  fi
}

ask_choose() {  # prompt, preset, options... -> echoes choice (default = first)
  local prompt="$1"; shift
  local preset="$1"; shift
  if [ -n "$preset" ]; then echo "$preset"; return 0; fi
  if [ "$INTERACTIVE" = "0" ]; then echo "$1"; return 0; fi
  if have gum; then
    printf '%s\n' "$@" | gum choose --header "$prompt"
  else
    local v
    PS3="$prompt: "
    select v in "$@"; do [ -n "$v" ] && { echo "$v"; break; }; done
  fi
}

ask_multi() {  # prompt, preset(csv), options... -> echoes csv of chosen
  local prompt="$1"; shift
  local preset="$1"; shift
  if [ -n "$preset" ]; then echo "$preset"; return 0; fi
  if [ "$INTERACTIVE" = "0" ]; then echo ""; return 0; fi
  [ "$#" -eq 0 ] && { echo ""; return 0; }
  if have gum; then
    printf '%s\n' "$@" | gum choose --no-limit --header "$prompt" | paste -sd, -
  else
    local v; read -r -p "$prompt (comma-separated, blank = none): " v; echo "$v"
  fi
}

confirm() {  # prompt -> exit status
  if [ "${BISR_ASSUME_YES:-0}" = "1" ]; then return 0; fi
  if have gum; then gum confirm "$1"; else
    local a; read -r -p "$1 [y/N]: " a; [[ "$a" == [yY]* ]]
  fi
}

# Front-end chooser: bash session (this script) or the Go TUI (tui/bisrde-tui).
# Skipped in non-interactive / --print-cmd mode so scripting + parity tests
# stay on the deterministic bash path.
TUI_BIN="$SCRIPT_DIR/tui/bisrde-tui"
choose_frontend() {
  [ "$PRINT_CMD" = "1" ] && return 0
  [ "$INTERACTIVE" = "0" ] && return 0
  local choice
  choice=$(ask_choose "Choose your interface" "${BISR_FRONTEND:-}" \
    "Bash interactive session" "Go TUI")
  if [ "$choice" = "Go TUI" ]; then
    if [ ! -x "$TUI_BIN" ]; then
      note "Go TUI not built yet."
      if have go && confirm "Build it now? (needs the Go toolchain)"; then
        ( cd "$SCRIPT_DIR/tui" && go build -o bisrde-tui . ) \
          || { note "Build failed — falling back to the bash session."; return 0; }
      else
        note "Build it with:  (cd tui && go build -o bisrde-tui .)"
        note "Falling back to the bash session."
        return 0
      fi
    fi
    cd "$SCRIPT_DIR"
    BISR_PROJECT_DIR="$SCRIPT_DIR" exec "$TUI_BIN"
  fi
  # "Bash interactive session" -> fall through to the bash flow below.
}

# ---------------------------------------------------------------------------
# Samplesheet introspection — extract groups / samples / contrast columns.
# Convention: column 1 = SampleID, column 2 = GroupID. Reserved metadata
# columns (Exclude, DisplayName) and the two ID columns are NOT contrasts.
# ---------------------------------------------------------------------------
ss_samples() { awk -F, 'NR>1 && NF>0 {print $1}' "$1" | sed '/^$/d'; }
ss_groups()  { awk -F, 'NR>1 && NF>0 {print $2}' "$1" | sed '/^$/d' | sort -u; }
ss_contrasts() {
  awk -F, 'NR==1{
    for (i=1;i<=NF;i++){
      h=$i; gsub(/\r/,"",h);
      if (h!="SampleID" && h!="GroupID" && h!="Exclude" && h!="DisplayName") print h
    }
  }' "$1"
}

# A samplesheet "has display names" when it carries a DisplayName column with
# at least one non-blank value. Exit 0 = yes, 1 = no (drives the #2b prompt).
ss_has_display() {
  awk -F, '
    NR==1{ for(i=1;i<=NF;i++){h=$i;gsub(/\r/,"",h);if(h=="DisplayName")c=i} if(!c)exit 1 }
    NR>1 && c { v=$c; gsub(/\r/,"",v); gsub(/^[ \t]+|[ \t]+$/,"",v); if(v!="")f=1 }
    END{ exit (f?0:1) }
  ' "$1"
}

# Write a working copy of the samplesheet with a DisplayName column filled from
# a group->label map, numbering replicates within each group ("<label> <n>").
# Args: <src> <group=label>...  Echoes the (absolute) working-copy path.
write_display_samplesheet() {
  local src="$1"; shift
  local mapf out pair
  mapf="$(mktemp)"; out="$(mktemp)"
  for pair in "$@"; do printf '%s\t%s\n' "${pair%%=*}" "${pair#*=}" >> "$mapf"; done
  awk -F, -v OFS=, -v mapf="$mapf" '
    BEGIN{ while((getline l < mapf)>0){ i=index(l,"\t"); lbl[substr(l,1,i-1)]=substr(l,i+1) } }
    { sub(/\r$/,"") }
    NR==1{ dc=0; for(i=1;i<=NF;i++) if($i=="DisplayName") dc=i
           print (dc? $0 : $0 ",DisplayName"); next }
    { g=$2; n[g]++; dn=((g in lbl)?lbl[g]:g) " " n[g]
      if(dc){ $dc=dn; print } else { print $0 "," dn } }
  ' "$src" > "$out"
  rm -f "$mapf"
  printf '%s\n' "$out"
}

# ===========================================================================
# Flow
# ===========================================================================
choose_frontend   # may exec into the Go TUI and never return
banner

counts=$(ask_input      "Counts TSV"        "assets/example_counts.tsv"        "${BISR_COUNTS:-}")
samplesheet=$(ask_input "Samplesheet CSV"   "assets/example_samplesheet.csv"   "${BISR_SAMPLESHEET:-}")
annotation=$(ask_choose "Genome annotation" "${BISR_ANNOTATION:-}" human mouse)
runid=$(ask_input       "Run ID"            "run_$(date +%Y%m%d_%H%M%S)"        "${BISR_RUNID:-}")
outdir=$(ask_input      "Output directory"  "./results"                        "${BISR_OUTDIR:-}")
brs=$(ask_input         "BRS ticket (blank = none)" ""                         "${BISR_BRS:-}")
idtype=$(ask_choose     "Gene ID type"      "${BISR_IDTYPE:-}" ensembl entrez symbol)
vlabels=$(ask_input     "Volcano labels — max genes to name (top N)" "10"           "${BISR_VOLCANO_LABELS:-}")

# v1.5.0 — interactive sample / group / contrast selection (only when the
# samplesheet is readable; otherwise these stay empty = no filtering).
excl_groups="" ; excl_samples="" ; incl_contrasts="" ; excl_contrasts=""
if [ -f "$samplesheet" ]; then
  note "Reading samplesheet for selection menus…"
  # Portable array fill (macOS ships bash 3.2, which lacks `mapfile`).
  _groups=();    while IFS= read -r _l; do _groups+=("$_l");    done < <(ss_groups "$samplesheet")
  _samples=();   while IFS= read -r _l; do _samples+=("$_l");   done < <(ss_samples "$samplesheet")
  _contrasts=(); while IFS= read -r _l; do _contrasts+=("$_l"); done < <(ss_contrasts "$samplesheet")

  excl_groups=$(ask_multi  "Exclude groups (optional)"  "${BISR_EXCLUDE_GROUPS:-}"  "${_groups[@]}")
  excl_samples=$(ask_multi "Exclude samples (optional)" "${BISR_EXCLUDE_SAMPLES:-}" "${_samples[@]}")

  con_mode=$(ask_choose "Contrasts to run" "${BISR_CONTRAST_MODE:-}" \
    "all" "include only a subset" "exclude a subset")
  case "$con_mode" in
    "include only a subset")
      incl_contrasts=$(ask_multi "Include these contrasts" "${BISR_INCLUDE_CONTRASTS:-}" "${_contrasts[@]}") ;;
    "exclude a subset")
      excl_contrasts=$(ask_multi "Exclude these contrasts" "${BISR_EXCLUDE_CONTRASTS:-}" "${_contrasts[@]}") ;;
  esac
  # Honor explicit contrast env vars even without a mode pick (parity with
  # the Go TUI's non-interactive behavior + convenient for scripting).
  [ -z "$incl_contrasts" ] && incl_contrasts="${BISR_INCLUDE_CONTRASTS:-}"
  [ -z "$excl_contrasts" ] && excl_contrasts="${BISR_EXCLUDE_CONTRASTS:-}"

  # v1.5.2 (#2b) — sample display names. If the samplesheet has no usable
  # DisplayName column, offer per-group labels (per-sample doesn't scale);
  # otherwise the R side auto-derives "<Group> <n>" with a warning. Prompt is
  # interactive-only: non-interactive/--print-cmd falls to the first option
  # (auto-derive), so the --samplesheet arg is unchanged (parity preserved).
  if ! ss_has_display "$samplesheet"; then
    dn_mode=$(ask_choose "Sample display names" "${BISR_DISPLAYNAME_MODE:-}" \
      "auto-derive from group + replicate (recommended)" \
      "enter a label per group")
    if [ "$dn_mode" = "enter a label per group" ]; then
      note "Tip: add a DisplayName column to your samplesheet for full control."
      _dn=()
      for _g in "${_groups[@]}"; do
        _dn+=("$_g=$(ask_input "Display label for group '$_g'" "$_g" "")")
      done
      samplesheet=$(write_display_samplesheet "$samplesheet" "${_dn[@]}")
      note "Using working samplesheet with display labels: $samplesheet"
    fi
  fi
fi

# ---------------------------------------------------------------------------
# Assemble the run_analysis.sh argument vector.
# ---------------------------------------------------------------------------
args=(--counts "$counts" --samplesheet "$samplesheet" --outdir "$outdir"
      --runid "$runid" --annotation "$annotation" --id-type "$idtype"
      --volcano-labels "$vlabels")
[ -n "$brs" ]            && args+=(--brs-ticket "$brs")
[ -n "$excl_samples" ]   && args+=(--exclude-samples "$excl_samples")
[ -n "$excl_groups" ]    && args+=(--exclude-groups "$excl_groups")
[ -n "$incl_contrasts" ] && args+=(--include-contrasts "$incl_contrasts")
[ -n "$excl_contrasts" ] && args+=(--exclude-contrasts "$excl_contrasts")

# ---------------------------------------------------------------------------
# Summary card.
# ---------------------------------------------------------------------------
summary_md=$(cat <<EOF
## Run configuration

| Setting        | Value |
|----------------|-------|
| Counts         | \`$counts\` |
| Samplesheet    | \`$samplesheet\` |
| Annotation     | $annotation |
| Run ID         | $runid |
| Output dir     | \`$outdir\` |
| BRS ticket     | ${brs:-—} |
| ID type        | $idtype |
| Volcano labels | $vlabels |
| Exclude groups | ${excl_groups:-—} |
| Exclude samples| ${excl_samples:-—} |
| Include contrasts | ${incl_contrasts:-—} |
| Exclude contrasts | ${excl_contrasts:-—} |
EOF
)

if have glow; then
  echo "$summary_md" | glow -
elif have gum; then
  gum style --border rounded --padding "1 2" --border-foreground "$BISR_ACCENT" "$summary_md"
else
  echo
  echo "$summary_md"
  echo
fi

cmd_preview="bash run_analysis.sh ${args[*]}"

if [ "$PRINT_CMD" = "1" ]; then
  echo "$cmd_preview"
  exit 0
fi

if ! confirm "Launch the pipeline with these settings?"; then
  note "Aborted — nothing was run."
  exit 0
fi

# ---------------------------------------------------------------------------
# Launch. We do NOT wrap in `gum spin` because the R pipeline streams rich
# cli progress output the user wants to watch live.
# ---------------------------------------------------------------------------
note "Launching… (full cli output follows)"
cd "$SCRIPT_DIR"

# Optional: screenshot the run-config card into the run's logs/ dir
# (charmbracelet freeze). Best-effort and gated — a no-op if freeze isn't
# installed. Runs after the cd so a relative --outdir resolves the same way
# run_analysis.sh resolves it. The final-summary panel is NOT screenshotted
# here because we exec the pipeline below (this process is replaced); it is
# preserved in <outdir>/logs/<ts>_session.log instead.
if have freeze && mkdir -p "${outdir%/}/logs" 2>/dev/null; then
  cfg_md="$(mktemp -t bisr_cfg.XXXXXX 2>/dev/null || echo /tmp/bisr_cfg.md)"
  printf '%s\n' "$summary_md" > "$cfg_md"
  freeze "$cfg_md" --language markdown -o "${outdir%/}/logs/config_card.png" \
    >/dev/null 2>&1 && note "Config card screenshot: ${outdir%/}/logs/config_card.png"
  rm -f "$cfg_md"
fi

exec bash run_analysis.sh "${args[@]}"
