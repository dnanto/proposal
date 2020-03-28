rule phy:
  input:
    root / "reg.fna"
  output:
    root / "phy" / "msa-1.fna",
    root / "phy" / "msa-1.log",
    root / "phy" / "run-1.treefile"
  params:
    config["cpu"],
    root / "phy" / "run-1"
  conda:
    "../envs/bio.yml"
  shell:
    """
      mafft --auto --adjustdirection --thread {params[0]:q} {input:q} > {output[0]:q} 2> {output[1]:q} && \
      rm -f {params[1]:q}.* && iqtree -s {output[0]:q} -pre {params[1]:q} -alrt 1000 -bb 1000 -bnni -nt {params[0]:q} > /dev/null;
    """

rule outlier:
  input:
    root / "phy" / "run-1.treefile",
    root / "phy" / "msa-1.fna"
  output:
    root / "phy" / "run-1.txt"
  params:
    config["alpha"],
    config["cpu"]
  conda:
    "../envs/R.yml"
  script:
    "../scripts/outliers.R"

rule rephy:
  input:
    root / "reg.fna",
    root / "phy" / "run-1.txt"
  output:
    root / "reg.fna.fai",
    root / "phy" / "msa-2.fna",
    root / "phy" / "msa-2.log",
    root / "phy" / "run-2.log"
  params:
    config["cpu"],
    root / "phy" / "run-2"
  shell:
    """
      rm -f {output[0]:q} && sed 's/_/|/g' {input[1]:q} | samtools faidx -r - {input[0]:q} | \
      mafft --auto --adjustdirection --thread {params[0]:q} - > {output[1]:q} 2> {output[2]:q} && \
      rm -f {params[1]:q}.* && iqtree -s {output[1]:q} -pre {params[1]:q} -alrt 1000 -bb 1000 -bnni -nt {params[0]:q} > /dev/null;
    """
