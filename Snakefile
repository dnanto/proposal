include: "rules/common.smk"

## targets ##

rule all:
    input:
        target_extract,
        target_phylo,
        target_beautify

## modules ##

include: "rules/targets.smk"
include: "rules/extract.smk"
include: "rules/phylo.smk"
include: "rules/beast.smk"
