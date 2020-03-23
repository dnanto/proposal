rule query:
  input:
    config["qry"]
  output:
    root / "qry.fna"
  shell:
    'cp "{input}" "{output}"'

rule local:
  input: 
    root / "qry.fna"
  output:
    root / "lcl.tsv"
  params:
    **config["align"],
    num_threads = config["cpu"]
  conda:
    "../envs/bio.yml"
  shell:    
    """
    blastn \
      -task {params.task:q} -db {params.db:q} -query {input:q} -out {output:q} -outfmt 7 \
      -perc_identity {params.perc_identity:q} -qcov_hsp_perc {params.qcov_hsp_perc:q} \
      -max_target_seqs {params.max_target_seqs:q} -num_threads {params.num_threads:q} \
      ;
    """

rule entry:
    input:
        root / "lcl.tsv"
    output:
        root / "lib.fna"
    params:
        config["align"]["db"]
    conda:
      "../envs/bio.yml"
    shell:
        """awk '/^[^#]/ {{ print $2; }}' {input:q} | uniq | blastdbcmd -db {params:q} -entry_batch - > {output:q};"""

rule glocal:
    input:
      config["qry"],
      root / "lib.fna"
    output:
      root / "glb.tsv"
    params:
      config["cpu"]
    conda:
      "../envs/bio.yml"
    shell:
      'glsearch36 -m 8CB -T {params[0]:q} {input[0]:q} {input[1]:q} > {output:q};'

rule feature:
  input:
    root / "glb.tsv"
  output:
    root / "src.tsv"
  params:
    **{**config["entrez"], **config["dates"]},
    qcov_identity = config["align"]["qcov_identity"]
  run:
    Entrez.email = params.email
    with open(input[0]) as file1, open(output[0], "w") as file2:
      key = "collection_date"
      keys = ("accessionversion", "title", "taxid", "collection_date")
      getter = itemgetter(*keys)
      print(*keys, sep = "\t", file = file2)
      subjects = (
          row["subject id"] for row in parse_outfmt7(file1) 
          if float(row["% identity"]) > params.qcov_identity
      )
      for batch in batchify(subjects, size=params.post_size):
          with Entrez.esummary(db=params.db, id=",".join(batch), retmode="json") as handle:
              for row in process_esummary(handle):
                  row[key] = normalize_date(row.get(key), params.formats)
                  print(*getter(row), sep = "\t", file = file2)

rule extract:
  input:
    root / "src.tsv",
    root / "glb.tsv",
    root / "lib.fna"
  output:
    target_extract
  run:
    # keep accessions with collection_date
    with open(input[0]) as file:
      getter = itemgetter("accessionversion", "collection_date", "taxid")
      reader = DictReader(file, delimiter="\t")
      meta = { ele[0]: ele[1:] for ele in map(getter, reader) if ele[1] != "NA" }
    data = {}
    # keep first entry for each accession
    with open(input[1]) as file:
      for row in parse_outfmt7(file):
        key = row["subject id"]
        if key in meta:
          data[key] = data.get(key, row)
    # index library
    idx = SeqIO.index_db(":memory:", input[2], "fasta")
    # output extract
    with open(output[0], "w") as file:
      for key, val in sorted(data.items(), key=lambda item: meta[item[0]][0]):
        start, end = int(val['s. start']), int(val['s. end'])
        rec = idx[key][start-1:end]
        rec.id = f"{val['subject id']}:{start}-{end}|{'|'.join(meta[key])}"
        rec.description = ""
        SeqIO.write(rec, file, "fasta")
