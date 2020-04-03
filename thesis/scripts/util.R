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
  separate(mut, c("REF", "ALT"), sep = " ") %>%
  filter(REF != ALT & REF != "N" & ALT != "N") %>%
  mutate(
    call = apply(., 1, function(row) {
      REF = row["REF"]
      ALT = row["ALT"]
      case_when(
        REF == "*" ~ "ins", ALT == "*" ~ "del",
        REF == "A" & ALT == "G" ~ "trs", REF == "G" & ALT == "A" ~ "trs",
        REF == "C" & ALT == "T" ~ "trs", REF == "T" & ALT == "C" ~ "trs",
        REF == "A" & ALT == "C" ~ "trv", REF == "C" & ALT == "A" ~ "trv",
        REF == "A" & ALT == "T" ~ "trv", REF == "T" & ALT == "A" ~ "trv",
        REF == "C" & ALT == "G" ~ "trv", REF == "G" & ALT == "C" ~ "trv",
        REF == "G" & ALT == "T" ~ "trv", REF == "T" & ALT == "G" ~ "trv",
        is_empty(intersect(code[[REF]], code[[ALT]])) ~ "dis",
        T ~ "sim"
      )
    })
  )

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
