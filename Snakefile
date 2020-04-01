include: "rules/common.smk"

## targets ##

rule all:
  input:
    root / "phylo" / "clock.str.rds",
    root / "phylo" / "clock.rlx.rds"

## modules ##

include: "rules/extract.smk"
include: "rules/phylo.smk"
