#!/usr/bin/env bash

# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# Before: >FN597015 Vitis vinifera, line PN40024, unoriented chromosome_1, chr1
# After: >chr1 FN597015|Vitis vinifera|PN40024|assembly12x|unoriented chromosome_1

outPrefix="Vvin-PN40024-12x-chr"

if [ ! -f VV_chr12x.fsa.gz ]; then
  echo -e "ERROR: can't find file VV_chr12x.fsa.gz\n" 1>&2
  exit 1
fi

zcat VV_chr12x.fsa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,","); gsub(" ","",b[4]); \
sub(" ","|",b[1]); split(b[2],c," "); sub(" ","",b[3]); \
printf ">"b[4]" "b[1]"|"c[2]"|assembly12x|"b[3]"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa
