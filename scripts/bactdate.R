#!/usr/bin/env RScript

hush <- function(expr) { suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr))) }

if (!pacman::p_isinstalled(BactDating)) pacman::p_install_gh("xavierdidelot/BactDating")

hush(library(lubridate))
hush(library(tidyverse))
hush(library(BactDating))

args <- commandArgs(trailingOnly = T)

nbIts <- as.integer(args[3])
path <- args[1]
prefix <- args[2]
tree <- loadCFML(prefix = tools::file_path_sans_ext(path))
tip.date <- decimal_date(ymd(str_extract(tree$tip.label, "\\d{4}-\\d{2}-\\d{2}")))

res1 <- bactdate(tree, tip.date, useRec = T, model = "poisson", nbIts = nbIts)
res2 <- bactdate(tree, tip.date, useRec = T, model = "relaxedgamma", nbIts = nbIts)

saveRDS(res1, str_c(prefix, "poi.rds", sep = "."))
saveRDS(res2, str_c(prefix, "rga.rds", sep = "."))
