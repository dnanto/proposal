#!/usr/bin/env python3

import json
import sys
from collections import OrderedDict
from datetime import datetime
from operator import itemgetter

from Bio import Entrez

params = dict(
    email="", 
    db="nuccore", 
    post_size=10, 
    formats = ("%d-%b-%Y", "%Y-%m-%d", "%b-%Y", "%Y"),
    to_fmt = "%Y-%m-%d",
    na_val = "NA",
    qcov_identity = 95
)

Entrez.email = params["email"]

def batchify(entries, size=10):
    batch = []
    for i, e in enumerate(entries, start=1):
        batch.append(e)
        if i % size == 0:
            yield batch
            batch = []

    if batch:
        yield batch

def process_esummary(file):
    obj = json.load(file)["result"]
    keys = ("accessionversion", "title", "taxid", "subtype", "subname")
    getter = itemgetter(*keys)
    for key in obj["uids"]:
        val = obj[key]
        row = OrderedDict((
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

def parse_outfmt7(file):
    fields = []
    for line in map(str.strip, file):
        if line.startswith("# Fields: "):
            fields = line[10:].split(", ")
        elif line and not line.startswith("#"):
            yield dict(zip(fields, line.split("\t")))


with sys.stdin as file1, sys.stdout as file2:
    key = "collection_date"
    args = (params["formats"], params["to_fmt"], params["na_val"])
    keys = ("accessionversion", "title", "taxid", "collection_date")
    getter = itemgetter(*keys)
    print(*keys, sep = "\t", file = file2)
    subjects = (
        row["subject id"] for row in parse_outfmt7(file1) 
        if float(row["% identity"]) > params["qcov_identity"]
    )
    for batch in batchify(subjects, size=params["post_size"]):
        print(batch)
        # with Entrez.esummary(db=params["db"], id=",".join(batch), retmode="json") as handle:
        #     for row in process_esummary(handle):
        #         row[key] = normalize_date(row.get(key), *args)
        #         print(*row.values(), sep = "\t", file = file2)
