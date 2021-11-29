## Required packages ##

packages <- c("ggplot2", "tidyverse", "data.table", "Rcpp", "rlang", "stringr", "cowplot", "RcppGSL", "qs", "viridis", "colorblindr", "readxl", "epiR", "crayon", "lubridate",
              "HDInterval","ggthemes", "ggpubr", "gridExtra", "tictoc", "lhs")

install.packages(setdiff(packages, rownames(installed.packages())))  

