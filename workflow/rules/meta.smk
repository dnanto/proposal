input_aln = root / ("cov.tsv" if config["mode"] else "glc.tsv")
input_seq = root / ("gen.fna" if config["mode"] else "lib.fna")

rule feature:
  input:
    input_aln
  output:
    root / "src.json"
  params:
    **config
  run:
    Entrez.email = params.email
    path = Path(input[0])
    field = "[Q]" if path.name == "cov.tsv" else "subject id"
    parser = parse_coords if path.name == "cov.tsv" else parse_outfmt7
    with path.open() as file1, open(output[0], "w") as file2:
      for batch in batchify(set(map(itemgetter(field), parser(file1))), size=params.post_size):
        with Entrez.esummary(db=params.edb, id=",".join(batch), retmode="json") as handle:
          file2.write(handle.read())

rule filter:
  input:
    input_aln,
    input_seq,
    root / "src.json"
  output:
    root / "ext.fna"
  params:
    config["qcov_idn_perc"],
    config["formats"]
  run:
    # file type
    path = Path(input[0])
    mode = path.name == "cov.tsv"
    field1 = "[Q]" if mode else "subject id"
    field2 = "[COV Q]" if mode else "% identity"
    parser = parse_coords if mode else parse_outfmt7
    # index library
    idx = SeqIO.index_db(":memory:", input[1], "fasta")
    # load metadata
    with open(input[2]) as file:
      meta = dict(chain.from_iterable(map(process_esummary, jsons(file))))
      date = { k: normalize_date(v.get("collection_date"), params[1]) for k, v in meta.items() }
    # filter records
    rec = {}
    with open(input[0]) as file:
      for row in parser(file):
        key = row[field1]
        if key not in rec and float(row[field2]) >= params[0] and date[key]:
          val = idx[key]
          if not mode:
            start, end = int(row["s. start"]), int(row["s. end"])
            val = val[start-1:end]
          val.description = val.description[len(val.id)+1:]
          val.id = f"{key}_{date[key]}_{meta[key]['taxid']}"
          rec[key] = val
    # output extract
    SeqIO.write(rec.values(), output[0], "fasta")
