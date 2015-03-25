#!/usr/bin/env bash

# Copyright (C) 2014-2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# Before: >FN597015 Vitis vinifera, line PN40024, unoriented chromosome_1, chr1
# After: >chr1 FN597015|Vitis vinifera|PN40024|assembly12x.0|unoriented chromosome_1

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_12x_v0_chroms_URGI"

if [ ! -f VV_chr12x.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_chr12x.fsa.gz\n" 1>&2
  exit 1
fi

date
echo "start ..."

zcat VV_chr12x.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,","); gsub(" ","",b[4]); \
sub(" ","|",b[1]); split(b[2],c," "); sub(" ","",b[3]); \
printf ">"b[4]" "b[1]"|"c[2]"|assembly12x.0|"b[3]"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
