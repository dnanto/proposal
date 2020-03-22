#!/usr/bin/env python3

import re
import sys
import time
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from datetime import datetime as dt
from pathlib import Path

from Bio import SeqIO
from bs4 import BeautifulSoup


def to_decimal_date(date):
    # https://stackoverflow.com/a/6451892
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year + 1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction

def taxa_tags(soup, path, dregex="(%Y-%m-%d)", dformat="%Y-%m-%d"):
    # taxa
    tag_tax = soup.new_tag("taxa", id="taxa")
    # alignment
    tag_aln = soup.new_tag("alignment", id="alignment", dataType="nucleotide")

    for rec in SeqIO.parse(path, "fasta"):
        # taxon
        value = to_decimal_date(dt.strptime(re.search(dregex, rec.id).group(1), dformat))
        tag_txn = soup.new_tag("taxon", id=rec.id)
        tag_txn.append(soup.new_tag("date", value=str(value), direction="forwards", units="years"))
        tag_tax.append(tag_txn)
        # sequence
        tag_seq = soup.new_tag("sequence")
        tag_seq.append(soup.new_tag("taxon", idref=rec.id))
        tag_seq.append(str(rec.seq).upper())
        tag_aln.append(tag_seq)
    
    return tag_tax, tag_aln

def gammaize(soup, model):
    soup.beast.siteModel.append(BeautifulSoup(f"""
        <gammaShape gammaCategories="4">
            <parameter id="siteModel_{model}.alpha" value="0.5" lower="0.0"/>
        </gammaShape>
    """, "lxml-xml"))
    soup.beast.operators.append(BeautifulSoup(f"""
        <scaleOperator scaleFactor="0.75" weight="0.1">
            <parameter idref="siteModel_{model}.alpha"/>
        </scaleOperator>
    """, "lxml-xml"))
    soup.beast.mcmc.joint.prior.append(BeautifulSoup(f"""
        <exponentialPrior mean="0.5" offset="0.0">
            <parameter idref="siteModel_{model}.alpha"/>
        </exponentialPrior>
    """, "lxml-xml"))
    soup.select_one("#fileLog").append(soup.new_tag("parameter", idref=f"siteModel_{model}.alpha"))


def invariantize(soup, model):
    soup.beast.siteModel.append(BeautifulSoup(f"""
        <proportionInvariant>
            <parameter id="siteModel_{model}.pInv" value="0.5" lower="0.0" upper="1.0"/>
        </proportionInvariant>
    """, "lxml-xml"))
    soup.beast.operators.append(BeautifulSoup(f"""
		<randomWalkOperator windowSize="0.75" weight="1" boundaryCondition="logit">
			<parameter idref="siteModel_{model}.pInv"/>
		</randomWalkOperator>
    """, "lxml-xml"))
    soup.beast.mcmc.joint.prior.append(BeautifulSoup(f"""
        <uniformPrior lower="0.0" upper="1.0">
            <parameter idref="siteModel_{model}.pInv"/>
        </uniformPrior>
    """, "lxml-xml"))
    soup.select_one("#fileLog").append(soup.new_tag("parameter", idref=f"siteModel_{model}.pInv"))

def parse_args(argv):
    parser = ArgumentParser(
        description="generate a BEAST input file based on an xml template...",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("path", type=FileType(), help="the path to the multiple sequence alignment FASTA file")
    parser.add_argument("tdir", help="the template directory", type=Path)
    parser.add_argument("model", help="the substitution model id")
    parser.add_argument("clock", help="the clock model")
    parser.add_argument("coalescent", help="the coalescent model")
    parser.add_argument("-dregex", help="the regular expression to extract the tip dates", default="(\d{4}-\d{2}-\d{2})")
    parser.add_argument("-dformat", help="the date format for tip dates", default="%Y-%m-%d")
    parser.add_argument("-gamma", action="store_true", help="the flag for site rate heterogeneity model")
    parser.add_argument("-invariant", action="store_true", help="the flag for invariant site model")
    return parser.parse_args(argv)

def main(argv):
    args = parse_args(argv[1:])
    
    with args.tdir.joinpath("model.xml").open() as file:
        model = BeautifulSoup(file, "xml").find("model", id=args.model)
        sub_model = model.select_one(f"#subModel_{args.model}")
        site_model = model.select_one("siteModel")
        operators = model.select_one("operators")
        prior = model.select_one("prior")
        log = model.select_one("log")
    
    # print(args.tdir.joinpath(f"{args.clock}-{args.coalescent}.xml"))
    with args.tdir.joinpath(f"{args.clock}-{args.coalescent}.xml").open() as file:
        soup = BeautifulSoup(file, "xml")
        tag_tax, tag_aln = taxa_tags(soup, args.path, args.dregex, args.dformat)
        soup.beast.insert(0, tag_tax)
        soup.beast.insert(1, tag_aln)
        soup.beast.insert(2, sub_model)
        soup.beast.insert(3, site_model)
        for ele in list(operators.children):
            soup.beast.operators.append(ele)
        for ele in list(prior.children):
            soup.beast.mcmc.joint.prior.append(ele)
        for ele in list(log.children):
            soup.select_one("#fileLog").append(ele)
        soup.treeDataLikelihood.partition.siteModel["idref"] += "_" + args.model
        if args.gamma:
            gammaize(soup, args.model)
        if args.invariant:
            invariantize(soup, args.model)

    print(soup)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
