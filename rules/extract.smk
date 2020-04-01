rule query:
  input:
    config["qry"]
  output:
    root / "qry.fna"
  run:
    copy2(input[0], output[0])

rule local:
  input:
    config["qry"]
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
      config["qry"],
      root / "lib.fna"
    output:
      root / "glc.tsv"
    params:
      config["cpu"]
    run:
      cmd = ("glsearch36", "-m", "8CB", "-T", str(params[0]), "-O", output[0], input[0], input[1])
      print(*cmd)
      run(cmd, stdout = DEVNULL)

rule feature:
  input:
    root / "glc.tsv"
  output:
    root / "src.json"
  params:
    **config
  run:
    Entrez.email = params.email
    with open(input[0]) as file1, open(output[0], "w") as file2:
      for batch in batchify(map(itemgetter("subject id"), parse_outfmt7(file1)), size=params.post_size):
        with Entrez.esummary(db=params.edb, id=",".join(batch), retmode="json") as handle:
          file2.write(handle.read())

rule extract:
  input:
    root / "lib.fna",
    root / "glc.tsv",
    root / "src.json"
  output:
    root / "reg.fna"
  params:
    config["qcov_idn_perc"],
    config["formats"]
  run:
    # index library
    idx = SeqIO.index_db(":memory:", input[0], "fasta")
    # load metadata
    with open(input[2]) as file:
      meta = dict(process_esummary(json.load(file)))
      date = { k: normalize_date(v.get("collection_date"), params[1]) for k, v in meta.items() }
    # filter records
    rec = {}
    with open(input[1]) as file:
      for row in parse_outfmt7(file):
        key = row["subject id"]
        if key not in rec and float(row["% identity"]) >= params[0] and date[key]:
          start, end = int(row["s. start"]), int(row["s. end"])
          val = idx[key][start-1:end]
          val.description = val.description[len(val.id)+1:]
          val.id = f"{key}_{date[key]}_{meta[key]['taxid']}"
          rec[key] = val
    # output extract
    SeqIO.write(rec.values(), output[0], "fasta")
