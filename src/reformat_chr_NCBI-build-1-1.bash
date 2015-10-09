#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html

# for assembled scaffolds:
# Before: >gi|225570623|ref|NC_012007.2|NC_012007 Vitis vinifera chromosome 1, reference assembly (based on 8x_WGS), whole genome shotgun sequence
# After: >chr1 NC_012007.2|Vitis vinifera|PN40024|assembly8x|chromosome_1

# for unassembled scaffolds:
# Before: >gi|221071967|ref|NW_002240929.1|VviUn_WGA355_1 Vitis vinifera genomic contig, reference assembly (based on 8x_WGS scaffold_3735)
# After: >chrUn_scaffold_3735 NW_002240929.1|Vitis vinifera|PN40024|assembly8x|chromosome_Un

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_8x_chroms_NCBI"

date
echo "start ..."

rm -f ${outPrefix}.fa
touch ${outPrefix}.fa
for chr in {1..19}; do
  if [ ! -f vvi_ref_chr${chr}.fa.gz ]; then
    echo -e "ERROR: can't find file vvi_ref_chr${chr}.fa.gz\n" 1>&2
    exit 1
  fi
  zcat vvi_ref_chr${chr}.fa.gz \
    | awk -v chr=${chr} 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,"|"); \
printf ">chr"chr" "b[4]"|Vitis vinifera|PN40024|assembly8x|chromosome_"chr"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
          >> ${outPrefix}.fa
done
chr="Un"
if [ ! -f vvi_ref_chr${chr}.fa.gz ]; then
  echo -e "ERROR: can't find file vvi_ref_chr${chr}.fa.gz\n" 1>&2
  exit 1
fi
zcat vvi_ref_chr${chr}.fa.gz \
  | awk -v chr=${chr} 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,"|"); \
split(b[5],c," "); gsub(")","",c[11]); \
printf ">chr"chr"_"c[11]" "b[4]"|Vitis vinifera|PN40024|assembly8x|chromosome_"chr"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
        >> ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
