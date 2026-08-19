#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
library_path <- Sys.getenv("RI_RLIB", "")
if (nzchar(library_path)) {
  dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(normalizePath(library_path), .libPaths()))
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = if (nzchar(library_path)) library_path else NULL)
}

required <- c(grf = "2.4.0", haven = "2.5.4")
for (package in names(required)) {
  installed <- requireNamespace(package, quietly = TRUE)
  correct <- installed && as.character(packageVersion(package)) == required[[package]]
  if (!correct) {
    remotes::install_version(
      package,
      version = required[[package]],
      lib = if (nzchar(library_path)) library_path else NULL,
      upgrade = "never"
    )
  }
}

for (package in names(required)) {
  observed <- as.character(packageVersion(package))
  if (observed != required[[package]]) {
    stop(sprintf("Expected %s %s, found %s", package, required[[package]], observed))
  }
  cat(sprintf("PASS: %s %s\n", package, observed))
}
