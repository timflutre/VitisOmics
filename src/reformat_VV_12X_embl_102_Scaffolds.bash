#!/usr/bin/env bash

# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# Before: >FN594950 Vitis vinifera, line PN40024, scaffold_0.assembly12x
# After: >scaffold_0 FN594950|Vitis vinifera|PN40024|assembly12x

outPrefix="Vvin-PN40024-12x-scaff"

if [ ! -f VV_12X_embl_102_Scaffolds.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_12X_embl_102_Scaffolds.fsa.gz\n" 1>&2
  exit 1
fi

zcat VV_12X_embl_102_Scaffolds.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,", "); split(b[3],c,"."); \
split(b[1],d," "); split(b[2],e," "); \
printf ">"c[1]" "d[1]"|"d[2]" "d[3]"|"e[2]"|"c[2]"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa
