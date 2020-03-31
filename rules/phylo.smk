rule msa:
  input:
    root / "reg.fna"
  output:
    root / "phylo" / "msa.fna",
    root / "phylo" / "msa.log",
  params:
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    """
      mafft --auto --adjustdirection --thread {params[0]:q} {input[0]:q} > {output[0]:q} 2> {output[1]:q};
    """

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

rule mcmc:
  input:
    root / "phylo" / "cfml.log"
  output:
    root / "phylo" / "clock.poi.rds",
    root / "phylo" / "clock.rga.rds"
  params:
    root / "phylo" / "clock",
    config["mcmc"]
  conda:
    "../envs/R.yml"
  shell:
    "./scripts/bactdate.R {input} {params[0]} {params[1]}"
