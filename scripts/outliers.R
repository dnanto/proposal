#!/usr/bin/env RScript

pacman::p_load(treedater)

path1 <- snakemake@input[[1]]
path2 <- snakemake@input[[2]]
alpha <- as.double(snakemake@params[[1]])
ncpu <- as.integer(snakemake@params[[2]])

tree <- read.tree(path1)
date <- sampleYearsFromLabels(tree$tip.label, regex = "\\d{4}-\\d{2}-\\d{2}")
slen <- ncol(read.dna(path2, format = "fasta", as.character = T))

x <- capture.output(tsim <- dater(tree, date, slen, clock = "strict", ncpu = ncpu))
x <- capture.output(tips <- outlierTips(tsim))

writeLines(
  as.character(dplyr::pull(dplyr::filter(tips, q >= alpha), taxon)),
  paste(tools::file_path_sans_ext(path1), "txt", sep = ".")
)
