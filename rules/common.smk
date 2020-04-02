import json
from datetime import datetime
from itertools import chain
from operator import itemgetter
from pathlib import Path
from shutil import copy2
from subprocess import DEVNULL, PIPE, Popen, run

from Bio import Entrez, SeqIO
from snakemake.utils import validate

## config ##

report: "../report/workflow.rst"
configfile: "conf/config.yml"
validate(config, "../conf/schema.yml")

## variables ##

root = Path(config["out"]) / Path(config["qry"]).stem
targets = [
    root / "phylo" / "clock.str.rds",
    root / "phylo" / "clock.rlx.rds"
]

## functions ##

def argify(conf, pfx="-"):
    arr = (ele.split("=", maxsplit = 1) for ele in conf)
    arr = ((pfx + ele[0], ele[1]) for ele in arr)
    return chain.from_iterable(arr)

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

def parse_coords(file):
    next(file)
    next(file)
    next(file)
    fields = (*next(file).split("\t"), "[R]", "[Q]")
    yield from (dict(zip(fields, line.strip().split("\t"))) for line in file)

def contextify(row):
    flank1, flank2 = int(row["q. start"]) - 1, int(row["query length"]) - int(row["q. end"])
    sstart, send, sstrand = int(row["s. start"]), int(row["s. end"]), row["subject strand"]
    sstart, send = (sstart - flank1, send + flank2) if sstrand == "plus" else (send - flank2, sstart + flank1)
    sstart = 1 if sstart < 1 else sstart
    return f'{row["subject acc.ver"]} {sstart}-{send} {sstrand}'

def normalize_date(val, formats, to_fmt = "%Y-%m-%d", na_val = None):
    result = na_val
    for fmt in formats:
        try:
            result = datetime.strptime(val, fmt).strftime(to_fmt)
        except:
            pass
        if result != na_val:
            break
    return result

def jsons(file):
    buffer = ""
    for line in file:
        buffer += line
        if line.startswith("}"):
            yield json.loads(buffer)
            buffer = ""

def process_esummary(obj):
    obj = obj["result"]
    for key in obj["uids"]:
        val = obj[key]
        meta = dict(zip(val["subtype"].split("|"), val["subname"].split("|")))
        yield val["accessionversion"], { **val, **meta }

