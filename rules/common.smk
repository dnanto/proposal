import json
from csv import DictReader
from datetime import datetime
from operator import itemgetter
from pathlib import Path

from Bio import Entrez, SeqIO
from snakemake.utils import validate


## setup ##

report: "../report/workflow.rst"
configfile: "conf/config.yml"
validate(config, "../conf/schema.yml")

## variables ##

root = Path(config["out"]) / Path(config["qry"]).stem

target = (
  root / "beast" / "rex-con.xml",
  root / "beast" / "rex-exp.xml",
  root / "beast" / "rln-con.xml",
  root / "beast" / "rln-exp.xml",
  root / "beast" / "str-con.xml",
  root / "beast" / "str-exp.xml"
)

## functions ##

def batchify(entries, size=10):
    batch = []
    for i, e in enumerate(entries, start=1):
        batch.append(e)
        if i % size == 0:
            yield batch
            batch = []

    if batch:
        yield batch

def parse_outfmt7(file):
    fields = []
    for line in map(str.strip, file):
        if line.startswith("# Fields: "):
            fields = line[10:].split(", ")
        elif line and not line.startswith("#"):
            yield dict(zip(fields, line.split("\t")))

def process_esummary(file):
    obj = json.load(file)["result"]
    keys = ("accessionversion", "title", "taxid", "subtype", "subname")
    getter = itemgetter(*keys)
    for key in obj["uids"]:
        val = obj[key]
        row = dict((
            *zip(keys, getter(val)),
            *zip(val["subtype"].split("|"), val["subname"].split("|"))
        ))
        del row["subtype"], row["subname"]
        yield row

def normalize_date(val, formats, to_fmt = "%Y-%m-%d", na_val = "NA"):
    result = na_val
    for fmt in formats:
        try:
            result = datetime.strptime(val, fmt).strftime(to_fmt)
        except:
            pass
        if result != na_val:
            break
    return result
