#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
import signal
from Bio import Entrez


def batchify(entries, size=10):
    batch = []
    for i, e in enumerate(entries, start=1):
        batch.append(e)
        if i % size == 0:
            yield batch
            batch = []

    if batch:
        yield batch


def parse_args(argv):
    parser = ArgumentParser(description="eutil", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "id", type=FileType(), help="the identifiers"
    )
    parser.add_argument(
        "-db", default="nuccore",
        help="the NCBI database"
    )
    parser.add_argument(
        "-eutil", default="efetch", help="the E-utility to use"
    )
    parser.add_argument(
        "-params", help="the space separated key=value pairs"
    )
    parser.add_argument(
        "-post-size", type=int, default=200,
        help="the number of records to post at a time"
    )
    parser.add_argument(
        "-email", default="",
        help="the e-mail to identify yourself to NCBI (for politeness reasons)"
    )
    args = parser.parse_args(argv)
    return args


def main(argv):
    args = parse_args(argv[1:])

    Entrez.email = args.email

    eutil = getattr(Entrez, args.eutil)
    params = dict(item.split("=") for item in args.params.split()) if args.params else {}

    with args.id as file:
        for batch in batchify(map(str.strip, file), size=args.post_size):
            with eutil(db=args.db, id=",".join(batch), **params) as handle:
                sys.stdout.write(handle.read())

    return 0


if __name__ == "__main__":
    signal.signal(signal.SIGINT, lambda x, y: None)
    sys.exit(main(sys.argv))
