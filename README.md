# dissertation

dissertation

# setup

```bash
# conda setup
# might need to uncomment the next line for some systems...
# export TMPDIR="$(pwd)/.tmp"
conda update -y -n base -c defaults conda
conda config --set env_prompt '({name}) '
conda install -y -c bioconda -c conda-forge snakemake-minimal

# example data & run...
mkdir -p data && efetch -id FJ643676.1 -db nuccore -format gb -mode text > data/FJ643676.1.gbk
./scripts/ffcds.py -split data/FJ643676.1.gbk && mv FJ643676.1/ data/
# remove --use-conda if you have the dependencies
snakemake -p --use-conda --cores all --config qry=data/FJ643676.1/ACU57030.fna --

# uncomment to clean if needed
# rm -rf "$TMPDIR"
```

NOTE: replace the config.yml parameters accordingly, such as the blast database
