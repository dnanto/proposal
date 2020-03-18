# dissertation

dissertation

# setup

```bash
# conda setup
# might need to uncomment the next line for some systems...
# export TMPDIR=~/tmp
conda update -n base -c defaults conda
conda env create --prefix ./env --file env.yml
conda activate ./env

# install babette: unavailable on conda, so gotta do this...
R --vanilla -e 'install.packages("babette", repos = "http://cran.r-project.org")'

# install BEAST2 packages...
export BEAST_PACKAGE_PATH="$(pwd)"/.pkg
bash -x ./sh/beastpkg.sh

# example data & run...
mkdir -p data && efetch -id FJ643676.1 -db nuccore -format gb -mode text > data/FJ643676.1.gbk
./py/ffcds.py -split data/FJ643676.1.gbk && mv FJ643676.1/ data/
snakemake -p --config qry=data/FJ643676.1/ACU57037.fna --
```

NOTE: replace the config.yml parameters accordingly, such as **bdb** (the blast database)
