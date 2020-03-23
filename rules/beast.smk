rule beautify:
  input:
    root / "msa.fna",
    root / "phy" / "run.log"
  output:
    target_beautify
  params:
    **config["beast"]
  conda:
    "../envs/bio.yml"
  shell:
    """
      line=( $(
        grep -v WARNING {input[1]:q} | \
          awk '{{$1=$1}}; f;/No./{{f=1}}' | \
          awk 'NR <= 132' | \
          grep -v "+R" | \
          sort -n -k 7 | \
          head -n 1 | \
          sed -e 's/+G4/+G/g' -e 's/+F//g'
      ) );
      echo ${{line[1]}};
      for ele in {output}; do
      model=( $(basename "$ele" | tr '[\-.]' ' ') );
      mkdir -p "$(dirname $ele)";
      ./scripts/beautify.py \
        -stem "${{ele/.xml/}}" \
        -len_mcmc {params.len_mcmc:q} -len_psss {params.len_psss:q} \
        -echo_mcmc {params.echo_mcmc:q} -echo_psss {params.echo_psss:q} \
        -path_steps {params.path_steps:q} \
        {input[0]:q} ./templates "${{line[1]}}" "${{model[0]}}" "${{model[1]}}" > $ele;
      done;
    """
