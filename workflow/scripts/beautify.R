#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(babette)


read_iqtree_model_test <- function(path)
{
  lines <- read_lines(path)
  lines <- lines[grep("^WARNING:", lines, invert = T)]
  idx <- grep("^ModelFinder", lines)
  read_table(lines[(idx+1):(idx+133)])
}

params <- as.list(snakemake@config)
params$root <- dirname(snakemake@input[[1]])

msa <- ape::read.dna(snakemake@input[[1]], "fasta")

best_model <-
  read_iqtree_model_test(snakemake@input[[2]]) %>%
  separate("Model", c("model", "mod1", "mod2"), sep = "\\+", extra = "drop", fill = "right") %>%
  filter(model %in% c("JC", "GTR", "HKY", "TN")) %>%
  mutate_at("model", recode, JC = "JC69", TrN = "TN93") %>%
  arrange(BIC) %>%
  pull(model) %>%
  first()

site_model <- create_site_model_from_name(best_model)
clock_models <- setNames(c(create_clock_model_strict, create_clock_model_rln), c("str", "rln"))
tree_priors <- setNames(c(create_tree_prior_ccp, create_tree_prior_cep), c("ccp", "cep"))

out <- file.path(params$root, "beast")
dir.create(out, showWarnings = F, recursive = T)

mcmc <- create_mcmc(
  chain_length = as.integer(params$chain_length),
  n_init_attempts = 10000
)

path_cal <- file.path(out, "cal.tsv")
write_tsv(
  data.frame(
    rownames(msa), 
    decimal_date(ymd(str_split(rownames(msa), "\\|", simplify = T)[,2]))
  ),
  path_cal, 
  col_names = F
)

paths <-
  expand.grid(names(clock_models), names(tree_priors)) %>%
  apply(1, function(row) {
    output <- file.path(out, str_c(row[1], row[2], sep = "-"), str_c("run", "xml", sep = "."))
    dir.create(dirname(output), showWarnings = F, recursive = T)
    create_beast2_input_file(
      file.path(params$root, "msa.fna"),
      output,
      site_model = site_model,
      clock_model = clock_models[[row[1]]](),
      tree_prior = tree_priors[[row[2]]](),
      mcmc = mcmc,
      tipdates_filename = path_cal
    )
    output
  })


run_mle <- str_c(
  '<run id="mcmc"',
  'spec="beast.gss.MultiThreadedNS"',
  str_c("threads", params$threads, sep = "="),
  str_c("chainLength", params$chain_length, sep = "="),
  str_c("subChainLength", params$sub_chain_length, sep = "="),
  str_c("particleCount", params$particle_count, sep = "="),
  'epsilon="1e-12">',
  sep = " "
)

for (path in paths)
{
  lines <- read_lines(path)
  # BEAST2 doesn't like this line and not sure why it's there...
  lines <- lines[grep('<taxonset idref="TaxonSet.G_VII_pre2003_msa"/>', lines, invert = T)]
  idx <- grep('<logger id="screenlog" logEvery="1000">', lines)
  lines[idx] <- str_replace(lines[idx], "1000", "100000")
  write_lines(lines, path)
  lines[grep("^<run", lines)] <- run_mle
  write_lines(lines, str_replace_all(path, "run", "mle"))
}
