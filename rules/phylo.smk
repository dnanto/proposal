rule msa:
  input:
    root / "reg.fna"
  output:
    root / "msa.fna",
    root / "msa.log",
  params:
    config["cpu"]
  conda:
    "../envs/phylo.yml"
  shell:
    """mafft --auto --adjustdirection --thread {params[0]:q} {input:q} > {output[0]:q} 2> {output[1]:q};"""

rule phy:
  input:
    root / "msa.fna"
  output:
    root / "phy" / "run.log"
  params:
    root / "phy" / "run",
    config["cpu"]
  conda:
    "../envs/phylo.yml"
  shell:
    """rm -f {params[0]:q}.* && iqtree -s {input:q} -pre {params[0]:q} -alrt 1000 -bb 1000 -bnni -nt {params[1]:q} > /dev/null;"""

rule beastify:
  input:
    root / "msa.fna",
    root / "phy" / "run.log"
  output:
    target
  params:
    **config["beast"]
  shell:
    """
      line=( $(
        grep -v WARNING {input[1]:q} | \
          awk '{{$1=$1}}; 261 >= NR && NR >= 130' | \
          grep -v "+R" | \
          sort -n -k 7 | \
          head -n 1 | \
          sed -e 's/+G4/+G/g' -e 's/+F//g'
      ) );
      echo ${{line[1]}};
      for ele in {output}; do
        model=( $(basename "$ele" | tr '[-.]' ' ') );
        mkdir -p "$(dirname $ele)";
        ./scripts/beautify.py \
          -stem "${{ele/.xml/}}" \
          -len_mcmc {params.len_mcmc:q} -len_psss {params.len_psss:q} \
          -echo_mcmc {params.echo_mcmc:q} -echo_psss {params.echo_psss:q} \
          -path_steps {params.path_steps:q} \
          {input[0]:q} ./templates "${{line[1]}}" ${{model[0]}} ${{model[1]}} > $ele;
      done;
    """
