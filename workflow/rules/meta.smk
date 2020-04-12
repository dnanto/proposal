rule esummary:
  input:
    root / f"ext-{mode}.fna"
  output:
    root / "meta.json"
  params:
    **config
  run:
    Entrez.email = params.email
    accs = (rec.id for rec in SeqIO.parse(input[0], "fasta"))
    with open(output[0], "w") as file:
      for batch in batchify(accs, size=params.post_size):
        with Entrez.esummary(db=params.edb, id=",".join(batch), retmode="json") as handle:
          file.write(handle.read())

rule ext_date:
  input:
    root / "meta.json"
  output:
    root / "date.tsv",
    root / "date.sed"
  params:
    config["formats"]
  run:
    with open(input[0]) as file1, open(output[0], "w") as file2, open(output[1], "w") as file3:
      meta = dict(chain.from_iterable(map(process_esummary, jsons(file1))))
      for key, val in meta.items():
        val = normalize_date(val.get("collection_date"), params[0])
        if val:
          print(key, val, sep="\t", file=file2)
          print(f"/^>/ s/>{key}/>{key}_{val}/", file=file3)
