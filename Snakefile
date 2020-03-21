include: "rules/common.smk"

## targets ##

rule all:
    input:
        target

## modules ##

include: "rules/extract.smk"
include: "rules/phylo.smk"
