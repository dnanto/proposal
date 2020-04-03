#!/usr/bin/env Rscript

hush <- function(expr) { suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr))) }

hush(library(lubridate))
hush(library(tidyverse))
hush(library(BactDating))

args <- commandArgs(trailingOnly = T)

path <- args[1]
model <- args[2]
iter <- as.integer(args[3])
thin <- as.integer(args[4])
out <- args[5]

tree <- loadCFML(prefix = tools::file_path_sans_ext(path))
tip.date <- decimal_date(ymd(str_extract(tree$tip.label, "\\d{4}-\\d{2}-\\d{2}")))

res <- bactdate(tree, tip.date, useRec = T, model = model, nbIts = iter, thin = thin, showProgress = T)

saveRDS(res, out)
