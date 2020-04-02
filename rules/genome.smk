rule genome:
  input:
    root / "ref.fna"
  output:
    root / "gen.fna"
  params:
    config["plen"],
    config["db"]
  run:
    n = len(SeqIO.read(input[0], "fasta"))
    low = n - n * params[0] / 100
    upp = n + n * params[0] / 100
    cmd = ("blastdbcmd", "-db", params[1], "-entry", "all", "-outfmt", "%a %l")
    print(cmd)
    with Popen(cmd, universal_newlines=True, stdout=PIPE) as proc:
      with proc.stdout as file:
        acc = [ele[0] for ele in map(str.split, file) if low <= float(ele[1]) <= upp]
    cmd = ("blastdbcmd", "-db", params[1], "-entry_batch", "-", "-out", output[0])
    with Popen(cmd, universal_newlines=True, stdin=PIPE) as proc:
      with proc.stdin as file:
        print(*acc, sep="\n", file=file)

rule nucmer:
  input:
    root / "ref.fna",
    root / "gen.fna"
  output:
    root / "aln.delta"
  params:
    config["cpu"],
    root / "aln"
  conda:
    "../envs/bio.yml"
  shell:
    "nucmer -t {params[0]:q} -p {params[1]:q} {input[0]:q} {input[1]:q}"

rule coords:
  input:
    root / "aln.delta"
  output:
    root / "cov.tsv"
  params:
    config["qcov_idn_perc"]
  conda:
    "../envs/bio.yml"
  shell:
    'show-coords -c -d -l -q -T {input[0]:q} > {output[0]:q}'
