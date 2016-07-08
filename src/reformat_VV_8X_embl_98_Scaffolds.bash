#!/usr/bin/env bash

# Copyright (C) 2016 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# Before: >CU459218 Vitis vinifera, line PN40024, chromosome chr18, scaffold_1
# After: >scaffold_1 CU459218|Vitis vinifera|PN40024|assembly8x|chr18

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_8x_scaffolds_EMBL_r98"

if [ ! -f VV_8X_embl_98_Scaffolds.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_8X_embl_98_Scaffolds.fsa.gz\n" 1>&2
  exit 1
fi

date
echo "start ..."

zcat VV_8X_embl_98_Scaffolds.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,", "); split(b[1],c," "); \
split(b[2],d," "); split(b[3],e," "); \
printf ">"b[4]" "c[1]"|"c[2]" "c[3]"|"d[2]"|assembly8x|"e[2]"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
