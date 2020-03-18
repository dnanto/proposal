from pathlib import Path
from csv import DictReader
from operator import itemgetter
from snakemake.utils import validate

# config

configfile: Path("conf") / "config.yml"
validate(config, Path("conf") / "schema.yml")

# setup

path = Path(config["qry"])
root = Path(config["out"]) / path.stem

def parse_outfmt7(file):
    fields = []
    for line in map(str.strip, file):
        if line.startswith("# Fields: "):
            fields = line[10:].split(", ")
        elif line and not line.startswith("#"):
            yield dict(zip(fields, line.split("\t")))

beautify_output = [
    root / "beast" / "cal.tsv",
    root / "beast" / "rln.ccp.mle.xml",
    root / "beast" / "rln.ccp.run.xml",
    root / "beast" / "rln.cep.mle.xml",
    root / "beast" / "rln.cep.run.xml",
    root / "beast" / "str.ccp.mle.xml",
    root / "beast" / "str.ccp.run.xml",
    root / "beast" / "str.cep.mle.xml",
    root / "beast" / "str.cep.run.xml"
]

# rules

rule all:
    input:
        beautify_output

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
        bdb=config["bdb"],
        idn=config["idn"],
        cov=config["cov"],
        thr=config["thr"]
    shell:    
        'blastn -task blastn -num_threads "{params.thr}" -max_target_seqs 20000 -perc_identity "{params.idn}" -qcov_hsp_perc "{params.cov}" -db "{params.bdb}" -query "{input}" -out "{output}" -outfmt 7;'

rule entry:
    input:
        root / "lcl.tsv"
    output:
        root / "lib.fna"
    params:
        config["bdb"]
    shell:
        """awk '/^[^#]/ {{ print $2; }}' {input} | uniq | blastdbcmd -db "{params}" -entry_batch - > "{output}";"""

rule glocal:
    input:
        config["qry"],
        root / "lib.fna"
    output:
        root / "glb.tsv"
    params:
        config["thr"]
    shell:
        'glsearch36 -m 8CB -T "{params[0]}" "{input[0]}" "{input[1]}" > "{output}";'

rule feature:
    input:
        root / "glb.tsv"
    output:
        root / "src.tsv",
        root / "src.json"
    params:
        pfx=root / "src",
        idn=config["idn"]
    shell:
        """awk -v idn="{params.idn}" '/^[^#]/ && $3 >= idn {{ print $2; }}' "{input}" | ./sh/feature.sh "{params.pfx}";"""

rule region:
    input:
        root / "src.tsv",
        root / "glb.tsv"
    output:
        root / "reg.txt",
        root / "reg.sed"
    run:
        # filter source
        with open(input[0]) as file:
            getter = itemgetter("accver", "collection_date", "taxid")
            reader = DictReader(file, delimiter="\t")
            meta = { ele[0]: ele[1:] for ele in map(getter, reader) if ele[1] != "NA" }
        # keep accessions with collection_date
        data = {}
        with open(input[1]) as file:
            for row in parse_outfmt7(file):
                key = row["subject id"]
                if key in meta:
                    data[key] = data.get(key, row)
        # output region and sed files
        with open(output[0], "w") as file1, open(output[1], "w") as file2:
            for key, val in sorted(data.items(), key=lambda item: meta[item[0]][0]):
                reg = f"{val['subject id']}:{val['s. start']}-{val['s. end']}"
                print(reg, file=file1)
                print(f"/^>/ s/{reg}/{reg}|{'|'.join(meta[key])}/g", file=file2)

rule extract:
    input:
        root / "lib.fna",
        root / "reg.txt",
        root / "reg.sed"
    output:
        root / "reg.fna"
    shell:
        """rm -f "{output}.fai" && samtools faidx "{input[0]}" -r "{input[1]}" | sed -E -f "{input[2]}" > "{output}";"""

rule msa:
    input:
        root / "reg.fna"
    output:
        root / "msa.fna",
        root / "msa.log",
    params:
        config["thr"]
    shell:
        """mafft --auto --adjustdirection --thread "{params[0]}" "{input}" > "{output[0]}" 2> "{output[1]}";"""

rule phyloml:
    input:
        root / "msa.fna"
    output:
        root / "phy" / "phy.log"
    params:
        pre=root / "phy" / "phy",
        thr=config["thr"]
    shell:
        """rm -f "{params.pre}."* && iqtree -s "{input}" -pre "{params.pre}" -alrt 1000 -bb 1000 -bnni -nt "{params.thr}" > /dev/null;"""

rule beautify:
    input:
        root / "msa.fna",
        root / "phy" / "phy.log"
    output:
        beautify_output
    script:
        "./R/beautify.R"
