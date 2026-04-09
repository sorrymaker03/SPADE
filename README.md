# SPADE

**S**patial **P**roximity **A**nalysis of **D**ifferential **E**xpression

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Overview

SPADE is an R package for spatial transcriptomics data analysis. It identifies cells located within specific distances from a reference cell type and performs differential expression analysis between cells in close proximity ("near") and those farther away ("far").

## Features

- Calculate distances between different cell types in spatial data
- Classify cells as "near", "mid", or "far" based on user-defined distance thresholds
- Perform differential expression analysis between near and far cells
- Filter out confounding signals using control cell populations
- Generate spatial visualizations of cell groups and distances
- Conduct GO enrichment and GSEA analysis on differentially expressed genes

## Installation

### From GitHub (recommended)

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SPADE
devtools::install_github("sorrymaker03/SPADE")

# Load the package
library(SPADE)
