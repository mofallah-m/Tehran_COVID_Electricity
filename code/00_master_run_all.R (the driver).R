## Set project root (edit as needed)
setwd("~/Desktop/TeIAS/Electricity")

## Fresh session friendliness
rm(list = ls()); graphics.off()

## Source in order
source("01_libraries_functions.R")
source("02_import_and_clean.R")
source("03_methods.R")
source("04_tests.R")
source("05_economic_cost.R")

message("All done. Figures and Tables have been written to ./Figures and ./Tables.")