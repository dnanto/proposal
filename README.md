# dissertation

dissertation

# setup

```bash
export TMPDIR=~/tmp # might need this for some systems...
mkdir -p .Rlib
R --vanilla -e 'install.packages(c("tidyverse", "lubridate"), lib=".Rlib", repos = "http://cran.r-project.org")'
conda env create --prefix ./env --file env.yml
conda activate ./env
mkdir -p data && efetch -id FJ643676.1 -db nuccore -format gb -mode text > data/FJ643676.1.gbk
./py/ffcds.py -split data/FJ643676.1.gbk && mv FJ643676.1/ data/
snakemake -p --config qry=data/FJ643676.1/ACU57037.fna --
```

note: replace the config.yml parameters accordingly, such as bdb (the blast database)
