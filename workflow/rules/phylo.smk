rule msa1:
  input:
    root / f"ext-{mode}.fna",
    root / "date.tsv",
    root / "date.sed"
  output:
    root / "msa-1.fna",
    root / "msa-1.log"
  params:
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    "rm -f {input[0]:q}.fai && "
    "cut -f 1 {input[1]:q} | "
    "samtools faidx -r - {input[0]:q} | "
    "sed -f {input[2]:q} | "
    "mafft --auto --adjustdirection --thread {params[0]:q} - > {output[0]:q} 2> {output[1]:q}"

rule recogub:
  input:
    root / "msa-1.fna"
  output:
    root / "gub.log",
    root / "gub.recombination_predictions.gff"
  params:
    root / "gub",
    config["reco_iter"],
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    "run_gubbins.py --prefix {params[0]:q} --iterations {params[1]:q} --threads {params[2]:q} {input[0]:q} > {output[0]:q}"

rule recobed:
  input:
    root / "msa-1.fna",
    root / "gub.recombination_predictions.gff",
  output:
    root / "rec.bed"
  run:
      msa = AlignIO.read(input[0], "fasta")

      with open(input[1]) as file:
        blocks = (ele for ele in file if not ele.startswith("#"))
        blocks = (list(map(int, ele.split("\t")[3:5])) for ele in blocks)
        blocks = sorted(blocks)
        blocks = ((ele[0] - 1, ele[1] + 1) for ele in contigify(blocks))
        blocks = chain.from_iterable(blocks)
        blocks = [0, *blocks, msa.get_alignment_length()]
        blocks = [[blocks[i], blocks[i+1]] for i in range(0, len(blocks), 2)]
        bsizes = ",".join(map(str, (ele[1] - ele[0] for ele in blocks)))
        bstart = ",".join((str(ele[0]) for ele in blocks))

      with open(output[0], "w") as file:
        for rec in msa:
            print(
                rec.id, 0, msa.get_alignment_length(), rec.id, 0, ".", 0, 0, 0,
                len(blocks), bsizes, bstart,
                sep="\t", file=file
            )

rule msa2:
  input:
    root / "msa-1.fna",
    root / "rec.bed"
  output:
    root / "msa-2.fna",
    root / "msa-2.log"
  params:
    config["cpu"]
  shell:
    "bedtools getfasta -fi {input[0]:q} -bed {input[1]:q} -fo - -name | "
    "mafft --auto --adjustdirection --thread {params[0]:q} - > {output[0]:q} 2> {output[1]:q}"

rule phy2:
  input:
    root / "msa-2.fna"
  output:
    root / "phy-2.treefile",
    root / "phy-2.log"
  params:
    root / "phy-2",
    config["cpu"]
  shell:
    "rm -f {params[0]:q}.* && "
    "iqtree -s {input[0]:q} -pre {params[0]:q} -alrt 1000 -bb 1000 -bnni -nt {params[1]:q} > /dev/null;"

rule mcmc_str:
  input:
    root / "gub.log"
  output:
    root / "clock.str.rds"
  params:
    config["mcmc_iter"],
    config["mcmc_thin"]
  conda:
    "../envs/R.yml"
  shell:
    "workflow/scripts/bactdate.R {input:q} poisson {params[0]:q} {params[1]:q} {output[0]:q}"

rule mcmc_rlx:
  input:
    root / "gub.log"
  output:
    root / "clock.rlx.rds"
  params:
    config["mcmc_iter"],
    config["mcmc_thin"]
  conda:
    "../envs/R.yml"
  shell:
    "workflow/scripts/bactdate.R {input:q} relaxedgamma {params[0]:q} {params[1]:q} {output[0]:q}"

rule beautify:
  input:
    root / "phy-2.log",
    root / "msa-2.fna"
  output:
    root / "rex-con.xml",
    root / "rex-exp.xml",
    root / "rln-con.xml",
    root / "rln-exp.xml",
    root / "str-con.xml",
    root / "str-exp.xml"
  params:
    **config
  conda:
    "../envs/bio.yml"
  shell:
    """
      line=( $(
        grep -v WARNING {input[0]:q} | awk '{{$1=$1;}} //1;' | \
          grep -A 131 '1 JC' | grep -v +R | sort -n -k 7 | head -n 1 | \
          sed -e 's/+G4/+G/g' -e 's/+F//g'
      ) );
      for ele in {output}; do
        model=( $(basename "$ele" | tr '[\-.]' ' ') );
        mkdir -p "$(dirname $ele)";
        workflow/scripts/beautify.py \
          -stem "${{model[0]}}-${{model[1]}}" \
          -len_mcmc {params.mcmc_len:q} -len_psss {params.psss_len:q} \
          -echo_mcmc {params.mcmc_echo:q} -echo_psss {params.psss_echo:q} \
          -path_steps {params.path_steps:q} \
          {input[1]:q} ./workflow/templates "${{line[1]}}" "${{model[0]}}" "${{model[1]}}" > $ele;
      done;
    """
