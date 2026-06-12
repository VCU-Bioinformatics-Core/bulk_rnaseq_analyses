# bisrDE/.Rprofile
# This subproject reuses the parent `differential_expression/` renv library
# to avoid duplicating Bioconductor + tidyverse deps.
#
# Strategy: prepend the parent's renv lib to .libPaths() WITHOUT bootstrapping
# a fresh sub-renv (which is what `source(parent_activate)` would do).
# All Imports listed in DESCRIPTION resolve via the parent lib.
local({
  parent_lib <- normalizePath(
    file.path(
      "..", "renv", "library", "R-4.2",
      paste0(R.version$arch, "-apple-darwin20")
    ),
    mustWork = FALSE
  )
  if (dir.exists(parent_lib)) {
    .libPaths(c(parent_lib, .libPaths()))
    message("bisrDE: parent renv library prepended to .libPaths()")
  } else {
    message(
      "bisrDE: parent renv library not found at ", parent_lib,
      "; falling back to system libs"
    )
  }
})
