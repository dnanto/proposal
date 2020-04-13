rule local:
  input:
    root / "ref-0.fna"
  output:
    root / "lcl.tsv"
  params:
    **config
  run:
    cmd = (
      "blastn", "-db", params.bdb, "-query", input[0], "-out", output[0], 
      "-outfmt", "7 std sstrand qlen", "-num_threads", str(params.cpu), *argify(params.blast)
    )
    taxidlist = params.get("taxidlist")
    cmd = (*cmd, "-taxidlist", taxidlist) if taxidlist else cmd
    print(*cmd)
    run(cmd) 

rule librarify:
  input:
    root / "lcl.tsv"
  output:
    root / "lib-0.fna"
  params:
    config["bdb"]
  run:
    cmd = ("blastdbcmd", "-db", params[0], "-entry_batch", "-")
    print(*cmd)
    with open(input[0]) as file1, open(output[0], "w") as file2:
      with Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdin as file:
          print(*map(contextify, parse_outfmt7(file1)), sep="\n", file=file)
        with proc.stdout as file:
          for rec in SeqIO.parse(file, "fasta"):
            rec.description = rec.description[len(rec.id)+1:]
            rec.id = rec.id.split(":", maxsplit=1)[0]
            SeqIO.write(rec, file2, "fasta")

rule glocal:
  input:
    root / "ref-0.fna",
    root / "lib-0.fna"
  output:
    root / "glc.tsv"
  params:
    config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:
    "glsearch36 -m 8CB -T {params[0]:q} {input[0]:q} {input[1]:q} > {output[0]:q}"

rule ext_out7:
  input:
    root / "glc.tsv",
    root / "lib-0.fna"
  output:
    root / "ext-0.fna"
  params:
    config["qcov_idn_perc"]
  conda:
    "../envs/bio.yml"
  shell:
    "awk -v cov={params[0]:q} '/^[^#]/ && $3 >= cov {{ print $2, $9, $10; }}' OFS='\t' {input[0]:q} | "
    "bedtools getfasta -fi {input[1]:q} -bed - -fo - | "
    "sed '/^>/ s/:[0-9]*-[0-9]*$//' > {output[0]:q}"
