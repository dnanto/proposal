#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)
conn <- if (args[1] == "-") file("stdin") else args[1]
col <- args[2]
orders <- if (length(args) >= 3) tail(args, -3) else c("dbY", "Ymd", "bY", "Y")

read_tsv(conn, col_types = cols(.default = "c")) %>%
  {
    if (col %in% names(.))
      mutate_at(., col, lubridate::parse_date_time, orders = orders, quiet = T) %>%
        mutate_at(col, strftime, format = "%Y-%m-%d")
    else .
  } %>%
  write.table("", row.names = F, quote = F, sep = "\t")
