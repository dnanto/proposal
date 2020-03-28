#!/usr/bin/env RScript

library(tidyverse)
library(treedater)

args <- commandArgs(trailingOnly = T)

path1 <- args[1]
path2 <- args[2]
alpha <- as.double(args[3])
ncpu <- as.integer(args[4])

tree <- read.tree(path1)
date <- sampleYearsFromLabels(tree$tip.label, regex = "\\d{4}-\\d{2}-\\d{2}")
slen <- read.dna(path2, format = "fasta", as.character = T) %>% ncol()

x <- capture.output(tsim <- dater(tree, date, slen, clock = "strict", ncpu = ncpu))
x <- capture.output(tips <- outlierTips(tsim) %>% filter(q <= 0.20) %>% pull(taxon))
cat(tips)
