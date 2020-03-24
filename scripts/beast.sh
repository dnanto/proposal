#!/usr/bin/env bash

burnin=${1:-10000000}
limit=${2:-0.5}

find out/*/beast -type f -name "*.xml" | \
while read -r ele; 
do 
    name="$(basename "$ele")";
    echo "time beast -overwrite "$ele" && time treeannotator -heights mean -burnin "$burnin" -limit "$limit" "${name/.xml/.trees}" "${name/.xml/.mcc.tree}"" | \
        qsub -V -d "$(dirname "$ele")" -N "${name//.xml}" -o "$root/${name/.xml/.out.log}" -e "$root/${name/.xml/.err.log}" -l nodes=1:ppn=32,walltime=48:00:00
done;
