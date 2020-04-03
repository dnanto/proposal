#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)
conn <- if (args[1] == "-") file("stdin") else args[1]
col1 <- args[2] # subtype
col2 <- args[3] # subname

read_tsv(conn, col_types = cols(.default = "c")) %>%
  bind_cols(
    .,
    bind_rows(
      apply(., 1, function(row) {
        key <- str_split(row[col1], "\\|")[[1]]
        val <- str_split(row[col2], "\\|")[[1]]
        data.frame(as.list(setNames(val, key)), stringsAsFactors = F)
      })
    )
  ) %>%
  select(-col1, -col2) %>%
  write.table("", row.names = F, quote = F, sep = "\t")
