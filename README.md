# dissertation

dissertation

# setup

```bash
conda env create --name dissertation -f env.yml
conda activate dissertation
mkdir -p data && efetch -id FJ643676.1 -db nuccore -format gb -mode text > data/FJ643676.1.gbk
./py/ffcds.py -split data/FJ643676.1.gbk && mv FJ643676.1/ data/
snakemake -p --config qry=data/FJ643676.1/ACU57037.fna --
```

note: replace the config.yml parameters accordingly, such as bdb (the blast database)
