package main

import "github.com/charmbracelet/lipgloss"

// Palette mirrors run_interactive.sh (pink primary, cyan accent).
var (
	colPrimary = lipgloss.Color("212")
	colAccent  = lipgloss.Color("86")
	colMuted   = lipgloss.Color("244")
	colOK      = lipgloss.Color("42")
	colErr     = lipgloss.Color("196")

	bannerStyle = lipgloss.NewStyle().
			Border(lipgloss.DoubleBorder()).
			BorderForeground(colPrimary).
			Foreground(colPrimary).
			Padding(1, 4).
			Margin(1, 0)

	paneStyle = lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(colAccent).
			Padding(1, 2).
			Margin(1, 0)

	mutedStyle = lipgloss.NewStyle().Foreground(colMuted)
	okStyle    = lipgloss.NewStyle().Foreground(colOK).Bold(true)
	errStyle   = lipgloss.NewStyle().Foreground(colErr).Bold(true)
	keyStyle   = lipgloss.NewStyle().Foreground(colAccent)
)
