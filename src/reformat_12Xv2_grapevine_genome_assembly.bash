#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# Before: >chr1 Assembly 12X.2
# After: >chr1 Vitis vinifera|PN40024|assembly12x.2

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_12x_v2_chroms_URGI"

if [ ! -f 12Xv2_grapevine_genome_assembly.fa.gz ]; then
  echo -e "ERROR: can't find file 12Xv2_grapevine_genome_assembly.fa.gz\n" 1>&2
  exit 1
fi

date
echo "start ..."

zcat 12Xv2_grapevine_genome_assembly.fa.gz \
  | awk 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b," "); \
printf ">"b[1]" Vitis vinifera|PN40024|assembly12x.2\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
  > ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
