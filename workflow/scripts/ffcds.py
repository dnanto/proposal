#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from pathlib import Path
from signal import SIG_DFL, SIGPIPE, signal

from Bio import SeqIO


def parse_cds(file):
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    protein_id = feature.qualifiers["protein_id"][0]
                    product = feature.qualifiers.get("product", ["n/a"])[0]
                    cds = feature.extract(record)
                    cds.id = protein_id
                    cds.description = f"{record.id}|{product}|{feature.location}"
                    yield protein_id, cds
                except:
                    pass


def parse_argv(argv):
    parser = ArgumentParser(
        description="extract CDS sequence records from a GenBank file",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("file", type=FileType(), help="the sequence file")
    parser.add_argument(
        "-split",
        action="store_true",
        help="the flag to split each cds into separate files",
    )
    args = parser.parse_args(argv)
    return args


def main(argv):
    args = parse_argv(argv[1:])

    with args.file as file:
        if args.split:
            root = Path(Path(file.name).stem)
            root.mkdir(exist_ok=True)
            for key, val in parse_cds(file):
                path = root.joinpath(key).with_suffix(".fna")
                print(path)
                SeqIO.write(val, path, "fasta")
        else:
            SeqIO.write((val for key, val in parse_cds(file)), sys.stdout, "fasta")

    return 0


if __name__ == "__main__":
    signal(SIGPIPE, SIG_DFL)
    sys.exit(main(sys.argv))
