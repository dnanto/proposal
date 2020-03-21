#!/usr/bin/env bash

prefix="$1"
entry="${2:-/dev/stdin}"

root="$(dirname "$0")"

{ 
  printf "accver\ttitle\ttaxid\tsubtype\tsubname\n";
  "$root"/../py/eutil.py "$entry" -eutil esummary -param retmode=json | \
    tee "$prefix.json" | \
    jq -r '.result | del(.uids) | map([.accessionversion, .title, .taxid, .subtype, .subname] | @tsv) | .[]';
} | \
  "$root"/../R/subtype.R - subtype subname 2> /dev/null | \
  "$root"/../R/lubridate.R - collection_date 2> /dev/null > "$prefix.tsv"
