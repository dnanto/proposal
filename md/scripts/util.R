code <- list(
  A = c("A"), C = c("C"), G = c("G"), T = c("T"),
  M = c("A", "C"), R = c("A", "G"), W = c("A", "T"), S = c("C", "G"), Y = c("C", "T"), K = c("G", "T"),
  V = c("A", "C", "G"), H = c("A", "C", "T"), D = c("A", "G", "T"), B = c("C", "G", "T"),
  N = c("A", "C", "G", "T")
)

muts <-
  outer(c(names(code), "*"), c(names(code), "*"), paste) %>%
  as.character() %>%
  enframe(name = NULL, value = "mut") %>%
  separate(mut, c("ref", "alt"), sep = " ") %>%
  filter(ref != alt & ref != "N" & alt != "N") %>%
  mutate(
    call = apply(., 1, function(row) {
      ref = row["ref"]
      alt = row["alt"]
      case_when(
        ref == "*" ~ "ins", alt == "*" ~ "del",
        ref == "A" & alt == "G" ~ "trs", ref == "G" & alt == "A" ~ "trs",
        ref == "C" & alt == "T" ~ "trs", ref == "T" & alt == "C" ~ "trs",
        ref == "A" & alt == "C" ~ "trv", ref == "C" & alt == "A" ~ "trv",
        ref == "A" & alt == "T" ~ "trv", ref == "T" & alt == "A" ~ "trv",
        ref == "C" & alt == "G" ~ "trv", ref == "G" & alt == "C" ~ "trv",
        ref == "G" & alt == "T" ~ "trv", ref == "T" & alt == "G" ~ "trv",
        is_empty(intersect(code[[ref]], code[[alt]])) ~ "dis",
        T ~ "sim"
      )
    })
  )

call_snp <- function(msa)
{
  ref <- msa[1, ]
  pos <- seq_along(ref)
  bind_rows(lapply(2:nrow(msa), function(idx) {
    alt <- msa[idx, ]
    pos <- pos[ref != '-' & alt != '-' & ref != alt]
    if (!is_empty(pos)) data.frame(idx = idx, pos = pos, ref = ref[pos], alt = alt[pos], stringsAsFactors = F)
  }))
}

call_ind <- function(msa)
{
  msa <- msa == "-"
  ref <- msa[1, ]
  bind_rows(lapply(2:nrow(msa), function(idx) {
    dff <- ref - msa[idx, ]
    bind_rows(
      str_locate_all(paste(abs(dff == +1), collapse = ""), "1+")[[1]] %>%
        as.data.frame() %>%
        mutate(idx = idx, call = "ins"),
      str_locate_all(paste(abs(dff == -1), collapse = ""), "1+")[[1]] %>%
        as.data.frame() %>%
        mutate(idx = idx, call = "del")
    )
  }))
}

is.int0 <- function(val) identical(val, integer(0))

tip.dates <- function(labs)
{
  decimal_date(ymd(str_extract(labs, "\\d{4}-\\d{2}-\\d{2}")))
}

read_snpsites <- function(path)
{
  lines <- read_lines(path)
  fields <- lines[startsWith(lines, "#CHROM")] %>% str_remove("^#")
  fields <- str_split(fields, "\t", simplify = T)[1,]
  df <- read_tsv(path, col_names = fields, col_types = cols(.default = "c"), comment = "#")
  alt <- str_split(df$ALT, ",", simplify = T)
  rownames(alt) <- df$POS
  pivot_longer(df, 10:ncol(df), names_to = "label", values_to = "index") %>%
    select(-CHROM, -ID, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    filter(index != "0") %>%
    mutate(ALT = apply(., 1, function(ele) alt[ele[["POS"]], as.integer(ele[["index"]])])) %>%
    mutate_at("POS", as.integer) %>%
    select(-index)
}

as_treedata <- function(res)
{
  d <- as.treedata(res1)
  treedata(phylo = d[[1]], data = as_tibble(d[[2]]))
}

plot_chronogram <- function(tree)
{
  tip.label <- tree@phylo$tip.label
  tip.date <- str_extract(tip.label, "\\d{4}-\\d{2}-\\d{2}")
  tip.taxid <- split(tip.label, str_extract(tip.label, "\\d+$"))
  tree@phylo <- groupOTU(tree@phylo, tip.taxid, "taxid")
  tree@phylo$tip.label <- str_remove(tip.label, "_\\d+$")
  ggtree(tree, mrsd = min(tip.date)) +
    aes(color = taxid) +
    geom_tiplab(linesize = 1, align = T, color = "black") +
    geom_range("height_0.95_HPD") +
    theme_tree2() +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_line(color="black", size = .25),
      panel.grid.minor.x = element_line(color="grey", size = .25)
    )
}

parse_mle_report <- function(path)
{
  lines <- read_lines(path)
  PS = "log marginal likelihood (using path sampling) from pathLikelihood.delta = "
  SS = "log marginal likelihood (using stepping stone sampling) from pathLikelihood.delta = "
  list(
    path = path,
    PS = as.double(str_sub(lines[startsWith(lines, PS)], nchar(PS) + 1)),
    SS = as.double(str_sub(lines[startsWith(lines, SS)], nchar(SS) + 1))
  )
}

parse_dot_attr <- function(val)
{
  tokens <-
    str_split(val, '(,)(?=(?:[^"]|"[^"]*")*$)') %>%
    lapply(str_split, '(=)(?=(?:[^"]|"[^"]*")*$)') %>%
    unlist() %>%
    str_trim() %>%
    str_remove('^"') %>%
    str_remove('"$')
  setNames(tokens[c(F,T)], tokens[c(T,F)])
}

read_dot <- function(path)
{
  lines <- read_lines(path) %>% str_trim()

  edge <- if(startsWith(lines[1], "digraph")) ">" else "-"

  e <-
    str_match(lines, str_c("(\\d+) -", edge, " (\\d+)")) %>%
    as.data.frame() %>%
    filter(complete.cases(.)) %>%
    select(-V1) %>%
    setNames(c("from", "to"))

  v <-
    str_match(lines, "(\\d+)\\[(.+)\\]") %>%
    as.data.frame() %>%
    filter(complete.cases(.)) %>%
    select(-V1) %>%
    setNames(c("id", "attr")) %>%
    bind_cols(map_df(lapply(.$attr, parse_attr), bind_rows))

  list(v, e)
}
