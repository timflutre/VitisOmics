#!/usr/bin/env bash

# Copyright (C) 2014-2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# Before: >FN594950 Vitis vinifera, line PN40024, scaffold_0.assembly12x
# After: >scaffold_0 FN594950|Vitis vinifera|PN40024|assembly12x.0

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_12x_v0_scaffolds_EMBL_r102"

if [ ! -f VV_12X_embl_102_Scaffolds.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_12X_embl_102_Scaffolds.fsa.gz\n" 1>&2
  exit 1
fi

date
echo "start ..."

zcat VV_12X_embl_102_Scaffolds.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,", "); split(b[3],c,"."); \
gsub(" ","",c[2]); split(b[1],d," "); split(b[2],e," "); \
printf ">"c[1]" "d[1]"|"d[2]" "d[3]"|"e[2]"|"c[2]".0\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
