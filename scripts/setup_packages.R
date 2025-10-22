############################################################
# Lightweight package setup helper
# - Installs missing packages via install.packages or pak if available
# - Optionally enforces minimum versions
############################################################

ensure_packages <- function(pkgs, min_versions = NULL) {
  # Set CRAN mirror if not already set
  if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cran.rstudio.com/"))
  }
  has_pak <- requireNamespace("pak", quietly = TRUE)
  for (i in seq_along(pkgs)) {
    pkg <- pkgs[i]
    need <- !requireNamespace(pkg, quietly = TRUE)
    if (!need && !is.null(min_versions) && !is.na(min_versions[i])) {
      cur <- as.character(utils::packageVersion(pkg))
      need <- utils::compareVersion(cur, min_versions[i]) < 0
    }
    if (need) {
      msg <- paste0("Installing package '", pkg, "'...")
      message(msg)
      if (has_pak) {
        pak::pak(pkg, ask = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
  }
}

# Example usage:
# ensure_packages(c("tidyverse","here","haven","nlme","imputeTS","doremi","lmerTest","sessioninfo"))

