#!/usr/bin/env bash

get_species_taxids.sh -t 10508 > 10508.txt
blastdbcmd -taxidlist 10508.txt -db ../nt/nt -outfmt "%a %T" > 10508.ssv
blastdbcmd -taxidlist 10508.txt -db ../nt/nt | \
  makeblastdb -in - \
    -dbtype nucl -title 10508 -out 10508 -logfile 10508.log \
    -taxid_map 10508.ssv -hash_index -parse_seqids -blastdb_version 5
