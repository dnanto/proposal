#!/usr/bin/env bash

burnin=${1:-10000000}
limit=${2:-0.5}

find out/*/beast -type f -name "*.xml" | \
while read -r ele; 
do 
    root="$(dirname $ele)";
    name="$(basename "$ele")";
    cd "$root" || exit
    echo "time beast -overwrite "$ele" && time treeannotator -heights mean -burnin "$burnin" -limit "$limit" "${name/.xml/.trees}" "${name/.xml/.mcc.tree}"" | \
        qsub -V -d "$root" -N "${name//.xml}" -o "${name/.xml/.out.log}" -e "${name/.xml/.err.log}" -l nodes=1:ppn=32,walltime=48:00:00
done;
