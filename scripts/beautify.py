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

def gammaize(soup):
    soup.beast.siteModel.append(BeautifulSoup("""
        <gammaShape gammaCategories="4">
            <parameter id="sitemodel.alpha" value="0.5" lower="0.0"/>
        </gammaShape>
    """, "lxml-xml"))
    soup.beast.operators.append(BeautifulSoup("""
        <scaleOperator scaleFactor="0.75" weight="0.1">
            <parameter idref="sitemodel.alpha"/>
        </scaleOperator>
    """, "lxml-xml"))
    soup.beast.mcmc.joint.prior.append(BeautifulSoup("""
        <exponentialPrior mean="0.5" offset="0.0">
            <parameter idref="sitemodel.alpha"/>
        </exponentialPrior>
    """, "lxml-xml"))
    soup.select_one("#fileLog").append(soup.new_tag("parameter", idref="sitemodel.alpha"))

def invariantize(soup):
    soup.beast.siteModel.append(BeautifulSoup("""
        <proportionInvariant>
            <parameter id="sitemodel.pInv" value="0.5" lower="0.0" upper="1.0"/>
        </proportionInvariant>
    """, "lxml-xml"))
    soup.beast.operators.append(BeautifulSoup("""
		<randomWalkOperator windowSize="0.75" weight="1" boundaryCondition="logit">
			<parameter idref="sitemodel.pInv"/>
		</randomWalkOperator>
    """, "lxml-xml"))
    soup.beast.mcmc.joint.prior.append(BeautifulSoup("""
        <uniformPrior lower="0.0" upper="1.0">
            <parameter idref="sitemodel.pInv"/>
        </uniformPrior>
    """, "lxml-xml"))
    soup.select_one("#fileLog").append(soup.new_tag("parameter", idref="sitemodel.pInv"))

def psss_tags(soup, path, **kwargs):
    with open(path) as file:
        for ele in list(BeautifulSoup(file, "xml").select_one("mle").children):
            soup.beast.append(ele)
        stem = kwargs["stem"]
        soup.select_one("marginalLikelihoodEstimator")["chainLength"] = kwargs["len_psss"]
        soup.select_one("marginalLikelihoodEstimator")["pathSteps"] = kwargs["path_steps"]
        soup.select_one("#MLELog")["logEvery"] = kwargs["echo_psss"]
        soup.select_one("#MLELog")["fileName"] = stem + ".mle.log"
        soup.select_one("pathSamplingAnalysis")["fileName"] = stem + ".mle.log"
        soup.select_one("pathSamplingAnalysis")["resultsFileName"] = stem + ".mle.result.log"
        soup.select_one("steppingStoneSamplingAnalysis")["fileName"] = stem + ".mle.log"
        soup.select_one("steppingStoneSamplingAnalysis")["resultsFileName"] = stem + ".mle.result.log"

def parse_args(argv):
    parser = ArgumentParser(
        description="generate a BEAST input file based on an xml template...",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    # taxa 
    parser.add_argument("msa", type=FileType(), help="the path to the multiple sequence alignment FASTA file")
    # model
    parser.add_argument("tdir", help="the template directory", type=Path)
    parser.add_argument("model", help="the substitution model id, add +I and +G for invariant and heterogeneity site models respectively")
    parser.add_argument("clock", help="the clock model")
    parser.add_argument("coalescent", help="the coalescent model")
    # parsing
    parser.add_argument("-dregex", help="the regular expression to extract the tip dates", default="(\d{4}-\d{2}-\d{2})")
    parser.add_argument("-dformat", help="the date format for tip dates", default="%Y-%m-%d")
    # mcmc
    parser.add_argument("-len_mcmc", type=int, help="the chain length for MCMC", default=10000000)
    parser.add_argument("-echo_mcmc", type=int, help="the sampling frequency for MCMC", default=10000)
    # psss
    parser.add_argument("-len_psss", type=int, help="the chain length for PS/SS", default=1000000)
    parser.add_argument("-echo_psss", type=int, help="the sampling frequency for PS/SS", default=1000)
    parser.add_argument("-path_steps", type=int, help="the number of path steps for PS/SS", default=100)
    # logging
    parser.add_argument("-stem", help="the output file stem", default="run")
    parser.add_argument("-echo", type=int, help="the echo frequency", default=0)

    return parser.parse_args(argv)

def main(argv):
    args = parse_args(argv[1:])

    with args.tdir.joinpath("model.xml").open() as file:
        model_id = args.model.split("+", maxsplit=1)[0]
        model = BeautifulSoup(file, "xml").find("model", id=model_id)
        sub_model = model.select_one("subModel")
        site_model = model.select_one("#sitemodel")
        operators = model.select_one("operators")
        prior = model.select_one("prior")
        log = model.select_one("log")
    
    with args.tdir.joinpath(f"{args.clock}-{args.coalescent}.xml").open() as file:
        soup = BeautifulSoup(file, "xml")
        # taxa
        tag_tax, tag_aln = taxa_tags(soup, args.msa, args.dregex, args.dformat)
        soup.beast.insert(0, tag_tax)
        soup.beast.insert(1, tag_aln)
        # model
        soup.beast.insert(2, sub_model)
        soup.beast.insert(3, site_model)
        if model.has_attr("operators"):
            for ele in list(operators.children):
                soup.beast.operators.append(ele)
        if model.has_attr("prior"):
            for ele in list(prior.children):
                soup.beast.mcmc.joint.prior.append(ele)
        if model.has_attr("log"):
            for ele in list(log.children):
                soup.select_one("#fileLog").append(ele)
        if "+G" in args.model:
            gammaize(soup)
        if "+I" in args.model:
            invariantize(soup)
        # MCMC
        soup.select_one("mcmc")["chainLength"] = args.len_mcmc
        soup.select_one("mcmc")["operatorAnalysis"] = args.stem + ".ops"
        soup.select_one("#fileLog")["logEvery"] = args.echo_mcmc
        soup.select_one("#fileLog")["fileName"] = args.stem + ".log"
        soup.select_one("logTree")["logEvery"] = args.echo_mcmc
        soup.select_one("logTree")["fileName"] = args.stem + ".trees"
        soup.select_one("#screenLog")["logEvery"] = args.echo
        # PS/SS
        psss_tags(soup, args.tdir.joinpath("psss.xml"), **vars(args))

    print(soup.prettify())

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
