# install_r_packages_for_visualization.R
# This script installs the required R packages for visualization of the results.

# Packages on CRAN
packages <- c("cowplot", "devtools", "ggseg", "ggsegGlasser", "glue", "rstatix", "see", "tidyverse")

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos="https://cloud.r-project.org/")
  }
}

sapply(packages, install_if_missing)

## LaCroixColoR from GitHub
 if (!requireNamespace("LaCroixColoR", quietly = TRUE)) {
    devtools::install_github("johannesbjork/LaCroixColoR")
}
