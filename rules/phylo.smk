rule msa:
  input:
    root / "ext.fna"
  output:
    root / "phylo" / "msa.fna"
  params:
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    """
      sed '/^>/ s/ .*//' {input[0]:q} | mafft --auto --adjustdirection --thread {params[0]:q} - > {output[0]:q} 2> {output[1]:q};
    """

rule snp:
  input:
    root / "phylo" / "msa.fna"
  output:
    root / "phylo" / "snp.vcf"
  shell:
    "snp-sites -v {input[0]:q} > {output[0]:q}"

rule phy:
  input:
    root / "phylo" / "msa.fna"
  output:
    root / "phylo" / "phy.treefile"
  params:
    root / "phylo" / "phy",
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    """
      rm -f {params[0]:q}.* && iqtree -s {input[0]:q} -pre {params[0]:q} -alrt 1000 -bb 1000 -bnni -nt {params[1]:q} > /dev/null;
    """

rule cfml:
  input:
    root / "phylo" / "phy.treefile",
    root / "phylo" / "msa.fna"
  output:
    root / "phylo" / "cfml.log"
  params:
    root / "phylo" / "cfml"
  conda:
    "../envs/bio.yml"
  shell:
    """ClonalFrameML {input[0]:q} {input[1]:q} {params[0]:q} -embranch true > {output[0]:q}"""

rule mcmc_str:
  input:
    root / "phylo" / "cfml.log"
  output:
    root / "phylo" / "clock.str.rds"
  params:
    config["iter"],
    config["thin"],
  conda:
    "../envs/R.yml"
  shell:
    './scripts/bactdate.R {input:q} "poisson" {params[0]:q} {params[1]:q} {output[0]:q}'

rule mcmc_rlx:
  input:
    root / "phylo" / "cfml.log"
  output:
    root / "phylo" / "clock.rlx.rds"
  params:
    config["iter"],
    config["thin"],
  conda:
    "../envs/R.yml"
  shell:
    './scripts/bactdate.R {input:q} "relaxedgamma" {params[0]:q} {params[1]:q} {output[0]:q}'
