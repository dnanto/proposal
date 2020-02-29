code <- list(
  A = c("A"), C = c("C"), G = c("G"), T = c("T"),
  M = c("A", "C"), R = c("A", "G"), W = c("A", "T"), S = c("C", "G"), Y = c("C", "T"), K = c("G", "T"),
  V = c("A", "C", "G"), H = c("A", "C", "T"), D = c("A", "G", "T"), B = c("C", "G", "T"),
  N = c("A", "C", "G", "T")
)

comp <- set_names(
  as.list(str_split("TAACGRYSWMKVHDBN", "", simplify = T)), 
  str_split("ATUGCYRSWKMBDHVN", "", simplify = T)
)

muts <-
  outer(c(names(code), "-"), c(names(code), "-"), paste) %>% 
  as.character() %>% 
  enframe(name = NULL, value = "mut") %>% 
  separate(mut, c("alt", "ref"), sep = " ", remove = F) %>% 
  mutate_at("mut", str_remove, " ") %>%
  mutate(
    call = apply(., 1, function(row) {
      ref = row["ref"]
      alt = row["alt"]
      case_when(
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
  ) %>%
  select(mut, call)

call_snp <- function(msa)
{
  unq <- apply(msa, 2, unique)
  idx <- which(sapply(unq, function(ele) !("-" %in% ele) && length(unique(ele)) != 1))
  tail(msa[ , idx], -1) %>% 
    as.data.frame() %>% 
    setNames(idx) %>% 
    rownames_to_column("id") %>% 
    pivot_longer(-id, names_to = "pos", names_ptypes = list(pos = integer()), values_to = "alt") %>%
    mutate(ref = msa[1, pos]) %>% 
    filter(ref != alt) %>%
    mutate(mut = str_c(ref, alt)) %>%
    merge(muts)
}

call_ind <- function(msa)
{
  apply(msa, 1, str_c, collapse = "") %>% 
    str_locate_all("-+") %>% 
    do.call(rbind, .) %>%
    as_tibble() %>%
    distinct() %>%
    apply(1, function(ele) {
      x <- msa[ , ele["start"]:ele["end"]]
      if (!is.null(ncol(x))) x <- apply(x, 1, str_c, collapse = "")
      c(ele["start"], x)
    }) %>%
    as.data.frame() %>%
    setNames(unlist(.[1, ])) %>%
    tail(-1) %>%
    rownames_to_column("id") %>%
    pivot_longer(-id, names_to = "pos", values_to = "alt") %>%
    merge(filter(., id == id[1]), by = c("pos"), suffixes = c("", ".y")) %>%
    filter(alt != alt.y) %>%
    select(pos, id, alt) %>%
    mutate(call = factor(if_else(str_detect(alt, "-"), "del", "ins"), c("ins", "del"))) %>%
    mutate(pos = as.integer(str_remove(pos, "\\..+")))
}
