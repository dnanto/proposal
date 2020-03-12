from pathlib import Path
from csv import DictReader
from operator import itemgetter
from collections import OrderedDict

path = Path(config["qry"])
root = Path(config["out"]) / Path(path.stem)

def parse_outfmt7(file):
	fields = []
	for line in map(str.strip, file):
		if line.startswith("# Fields: "):
			fields = line[10:].split(", ")
		elif line and not line.startswith("#"):
			yield OrderedDict(zip(fields, line.split("\t")))

rule all:
    input:
        root / "phy.treefile"

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
        config["db"]
    shell:    
        'blastn -task blastn -num_threads 8 -max_target_seqs 20000 -db "{params[0]}" -query "{input}" -out "{output}" -outfmt 7;'

rule entry:
    input:
        root / "lcl.tsv"
    output:
        root / "lib.fna"
    params:
        config["db"]
    shell:
        """awk '/^[^#]/ {{ print $2; }}' {input} | uniq | blastdbcmd -db "{params}" -entry_batch - > "{output}";"""

rule glocal:
    input:
        config["qry"],
        root / "lib.fna"
    output:
        root / "glb.tsv"
    shell:
        'glsearch36 -m 8CB -T 8 "{input[0]}" "{input[1]}" > "{output}";'

rule feature:
    input:
        root / "glb.tsv"
    output:
        root / "src.tsv",
        root / "src.json"
    params:
        pfx=root / "src",
        pct=config["pct"]
    shell:
        """awk -v pct="{params.pct}" '/^[^#]/ && $3 >= $pct {{ print $2; }}' "{input}" | ./sh/feature.sh "{params.pfx}";"""

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
            meta = { ele[0]: "|".join(ele[1:]) for ele in map(getter, reader) if ele[1] != "NA" }
        # keep accessions with collection_date
        data = {}
        with open(input[1]) as file:
            for row in parse_outfmt7(file):
                key = row["subject id"]
                if key in meta:
                    data[key] = data.get(key, row)
        # output region and sed files
        with open(output[0], "w") as file1, open(output[1], "w") as file2:
            for key, val in data.items():
                reg = f"{val['subject id']}:{val['s. start']}-{val['s. end']}"
                print(reg, file=file1)
                print(f"/^>/ s/{reg}/{reg}|{meta[key]}/g", file=file2)

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
        root / "msa.log"
    shell:
        """mafft --auto --adjustdirection --thread -1 "{input}" > "{output[0]}" 2> "{output[1]}";"""

rule phyloml:
    input:
        root / "msa.fna"
    output:
        root / "phy.treefile"
    params:
        root / "phy"
    shell:
        """iqtree -s "{input}" -pre "{params}" -alrt 1000 -bb 1000 -bnni -nt AUTO > /dev/null;"""
