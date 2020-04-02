rule local:
  input:
    root / "cds.fna"
  output:
    root / "lcl.tsv"
  params:
    **config
  run:
    cmd = (
      "blastn", "-db", params.db, "-query", input[0], "-out", output[0], 
      "-outfmt", "7 std sstrand qlen", "-num_threads", str(params.cpu), *argify(params.blast)
    )
    print(*cmd)
    run(cmd) 

rule entry:
  input:
    root / "lcl.tsv"
  output:
    root / "lib.fna"
  params:
    config["db"]
  run:
    cmd = ("blastdbcmd", "-db", params[0], "-entry_batch", "-")
    print(*cmd)
    with open(input[0]) as file1, open(output[0], "w") as file2:
      with Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdin as file:
          print(*map(contextify, parse_outfmt7(file1)), sep="\n", file = file)
        with proc.stdout as file:
          for rec in SeqIO.parse(file, "fasta"):
            rec.description = rec.description[len(rec.id)+1:]
            rec.id = rec.id.split(":", maxsplit = 1)[0]
            SeqIO.write(rec, file2, "fasta")

rule glocal:
    input:
      root / "cds.fna",
      root / "lib.fna"
    output:
      root / "glc.tsv"
    params:
      config["cpu"]
    run:
      cmd = ("glsearch36", "-m", "8CB", "-T", str(params[0]), "-O", output[0], input[0], input[1])
      print(*cmd)
      run(cmd, stdout = DEVNULL)
