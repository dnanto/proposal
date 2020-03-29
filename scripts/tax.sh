#!/usr/bin/env bash

root="data/tax"
mkdir -p "$root"
cd "$root" || exit

rsync -chavP --stats "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" .
tar -xf taxdump.tar.gz names.dmp

rm -f tax.sdb
sqlite3 tax.sdb << EOF
CREATE TABLE "name"(
  "tax_id" INTEGER NOT NULL,
  "name_txt" TEXT NOT NULL,
  "unique_name" TEXT,
  "name_class" TEXT NOT NULL
);
EOF

sed -e $'s/\t\|//g' names.dmp | sqlite3 -separator $'\t' tax.sdb ".import /dev/stdin name"
sqlite3 tax.sdb 'CREATE INDEX "name.index" ON "name" ("tax_id");'

rm names.dmp
