# bisrDE terminal demo — regenerate with `just vhs`.
#
# VHS (github.com/charmbracelet/vhs) records this SCRIPT, so the demo GIF is
# reproducible: every release re-records identically instead of someone
# hand-performing a run. Same Charm ecosystem as the gum/glow/freeze the
# launcher already uses.
#
# Requires: vhs (brew install vhs) and the bundled example data.

Output ../docs/demo.gif

Set Shell bash
Set FontSize 15
Set Width 1300
Set Height 760
Set Padding 18
Set Theme "Catppuccin Mocha"
Set TypingSpeed 45ms

# --- the one discoverable surface -------------------------------------------
Type "just --list"
Enter
Sleep 4s

# --- assemble a run without executing it (parity contract) ------------------
Type "clear"
Enter
Sleep 500ms
Type "just print-cmd"
Enter
Sleep 5s

# --- a real end-to-end run on the bundled example data ----------------------
Type "clear"
Enter
Sleep 500ms
Type "just demo"
Enter
# Long enough to capture the live progress bar advancing through the
# per-comparison sub-steps (DESeq2 -> annotate -> volcano -> heatmaps -> GSEA).
Sleep 100s

# --- what the run produced --------------------------------------------------
Type "clear"
Enter
Sleep 500ms
Type "ls results/"
Enter
Sleep 3s
