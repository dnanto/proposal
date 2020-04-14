rule subset:
  input:
    root / "ref-1.fna"
  output:
    root / "lib-1.fna"
  params:
    plen = config["plen"],
    bdb = config["bdb"],
    taxidlist = config.get("taxidlist")
  run:
    n = len(SeqIO.read(input[0], "fasta"))
    low = n - n * params.plen / 100
    upp = n + n * params.plen / 100
    args = ("-taxidlist", params.taxidlist) if params.taxidlist else ("-entry", "all")
    cmd = ("blastdbcmd", "-db", params.bdb, "-outfmt", "%a %l", *args)
    with Popen(cmd, universal_newlines=True, stdout=PIPE) as proc:
      with proc.stdout as file:
        acc = [ele[0] for ele in map(str.split, file) if low <= float(ele[1]) <= upp]
    cmd = ("blastdbcmd", "-db", params.bdb, "-entry_batch", "-", "-out", output[0])
    with Popen(cmd, universal_newlines=True, stdin=PIPE) as proc:
      with proc.stdin as file:
        print(*acc, sep="\n", file=file)

rule nucmer:
  input:
    root / "ref-1.fna",
    root / "lib-1.fna"
  output:
    root / "coor.delta"
  threads:
    8
  params:
    root / "coor"
  conda:
    "../envs/bio.yml"
  shell:
    "nucmer -t {threads} -p {params[0]:q} {input[0]:q} {input[1]:q}"

rule coords:
  input:
    root / "coor.delta"
  output:
    root / "coor.tsv"
  conda:
    "../envs/bio.yml"
  shell:
    "show-coords -c -d -l -q -T {input[0]:q} > {output[0]:q}"

rule ext_coor:
  input:
    root / "coor.tsv",
    root / "lib-1.fna"
  output:
    root / "ext-1.fna"
  params:
    config["qcov_idn_perc"]
  conda:
    "../envs/bio.yml"
  shell:
    "awk -v cov={params[0]:q} 'NR >= 5 && $11 >= cov {{ print $15, $3, $4 }}' OFS='\t' {input[0]:q} | "
    "bedtools getfasta -fi {input[1]:q} -bed - -fo - | "
    "sed '/^>/ s/:[0-9]*-[0-9]*$//' > {output[0]:q}"
