include: "rules/common.smk"

## rule targets ##

rule all:
  input:
    targets

## rule init ##

rule query:
  input:
    config["qry"]
  output:
    root / ("ref.fna" if config["mode"] else "cds.fna")
  run:
    copy2(input[0], output[0])

## rule modules ##

include: "rules/gene.smk"
include: "rules/genome.smk"
include: "rules/meta.smk"
include: "rules/phylo.smk"
