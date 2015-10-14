#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html

# Localized sequence:
# Before: >gi|221070425|ref|NW_002238163.1|Vvi1_WGA9_1 Vitis vinifera chromosome 1 genomic contig, reference assembly (based on 8x_WGS scaffold_75)
# After: >scaffold_75 NW_002238163.1|Vitis vinifera|PN40024|assembly8x|chr1

# Unplaced sequence:
# Before: >gi|221071967|ref|NW_002240929.1|VviUn_WGA355_1 Vitis vinifera genomic contig, reference assembly (based on 8x_WGS scaffold_3735)
# After: >scaffold_3735 NW_002240929.1|Vitis vinifera|PN40024|assembly8x|unplaced chrUn

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_8x_scaffolds_NCBI"

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
split($0,a,"\n"); split(a[1],b,"|"); split(b[5],c," "); \
gsub(")","",c[13]); \
printf ">"c[13]" "b[4]"|Vitis vinifera|PN40024|assembly8x|chr"chr"\n";
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
split($0,a,"\n"); split(a[1],b,"|"); split(b[5],c," "); \
gsub(")","",c[11]); \
printf ">"c[11]" "b[4]"|Vitis vinifera|PN40024|assembly8x|chr"chr"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
        >> ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
