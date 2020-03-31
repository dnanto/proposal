include: "rules/common.smk"

## targets ##

rule all:
  input:
    root / "phylo" / "clock.poi.rds",
    root / "phylo" / "clock.rga.rds"

## modules ##

include: "rules/extract.smk"
include: "rules/phylo.smk"
