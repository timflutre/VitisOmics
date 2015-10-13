#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html

# Localized sequence:
# Before: >gi|357741471|ref|NW_003723981.1| Vitis vinifera cultivar PN40024 chromosome 1 genomic scaffold, 12X scaffold_11.assembly12x, whole genome shotgun sequence
# After: >scaffold_11 NW_003723981.1|Vitis vinifera|PN40024|assembly12x.0|localized chr1

# Unlocalized sequence:
# Before: >gi|357737298|ref|NW_003724154.1| Vitis vinifera cultivar PN40024 chromosome 1 unlocalized genomic scaffold, 12X scaffold_155.assembly12x, whole genome shotgun sequence
# After: >scaffold_155 NW_003724154.1|Vitis vinifera|PN40024|assembly12x.0|unlocalized chr1

# Unplaced sequence:
# Before: >gi|357737221|ref|NW_003724231.1| Vitis vinifera cultivar PN40024 unplaced genomic scaffold, 12X scaffold_226.assembly12x, whole genome shotgun sequence
# After: >scaffold_266 NW_003724231.1|Vitis vinifera|PN40024|assembly12x.0|unplaced chrUn

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

outPrefix="VITVI_PN40024_12x_v0_scaffolds_NCBI"

date
echo "start ..."

rm -f ${outPrefix}.fa
touch ${outPrefix}.fa

for chr in {1..19}; do
  if [ ! -f vvi_ref_12X_chr${chr}.fa.gz ]; then
    echo -e "ERROR: can't find file vvi_ref_12X_chr${chr}.fa.gz\n" 1>&2
    exit 1
  fi
  zcat vvi_ref_12X_chr${chr}.fa.gz \
    | awk -v chr=${chr} 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,"|"); split(b[5],c," "); \
if(match(b[5],"unlocalized")==0){split(c[10],d,"\\."); \
printf ">"d[1]" "b[4]"|Vitis vinifera|PN40024|assembly12x.0|localized chr"chr"\n"} \
else{split(c[11],d,"\\."); \
printf ">"d[1]" "b[4]"|Vitis vinifera|PN40024|assembly12x.0|unlocalized chr"chr"\n"}; \
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
          >> ${outPrefix}.fa
done

chr="Un"
if [ ! -f vvi_ref_12X_chr${chr}.fa.gz ]; then
  echo -e "ERROR: can't find file vvi_ref_12X_chr${chr}.fa.gz\n" 1>&2
  exit 1
fi
zcat vvi_ref_12X_chr${chr}.fa.gz \
  | awk -v chr=${chr} 'BEGIN{RS=">"} {if(NF==0)next; \
split($0,a,"\n"); split(a[1],b,"|"); \
split(b[5],c," "); split(c[9],d,"\\."); \
printf ">"d[1]" "b[4]"|Vitis vinifera|PN40024|assembly12x.0|unplaced chr"chr"\n";
for(i=2;i<length(a);++i)printf a[i]"\n"; \
printf a[length(a)]}' \
        >> ${outPrefix}.fa

rm -f ${outPrefix}.fa.gz
gzip ${outPrefix}.fa

echo "done!"
date
