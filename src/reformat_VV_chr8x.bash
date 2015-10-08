#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# Before: >CU462738 Vitis vinifera, line PN40024, chromosome_1, chr1
# After: >chr1 CU462738|Vitis vinifera|PN40024|assembly8x|chromosome_1

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_8x_chroms_URGI"

if [ ! -f VV_chr8x.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_chr8x.fsa.gz\n" 1>&2
  exit 1
fi

date
echo "start ..."

zcat VV_chr8x.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,","); split(b[1],c," "); \
if(length(c)==1){printf ">"c[1]"\n"} \
else{gsub(" ","",b[4]); gsub(" ","",b[3]); \
printf ">"b[4]" "c[1]"|Vitis vinifera|PN40024|assembly8x|"b[3]"\n"};
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
