rule msa:
  input:
    root / "reg.fna"
  output:
    root / "msa.fna",
    root / "msa.log",
  params:
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    """mafft --auto --adjustdirection --thread {params[0]:q} {input:q} > {output[0]:q} 2> {output[1]:q};"""

rule phy:
  input:
    root / "msa.fna"
  output:
    target_phylo
  params:
    root / "phy" / "run"
  conda:
    "../envs/bio.yml"
  shell:
    """rm -f {params[0]:q}.* && iqtree -s {input:q} -pre {params[0]:q} -alrt 1000 -bb 1000 -bnni -nt AUTO > /dev/null;"""
