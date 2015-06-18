#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Gautier SARAH [cre,aut]

# http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

date
echo "start ..."

if [ ! -f VITVI_PN40024_12x_v2_chroms_URGI.fa.gz ]; then
  echo -e "ERROR: missing file VITVI_PN40024_12x_v2_chroms_URGI.fa.gz\n" 1>&2
  exit 1
fi
zcat VITVI_PN40024_12x_v2_chroms_URGI.fa.gz >VITVI_PN40024_12x_v2_chroms_URGI.tmp.fa
bowtie2-build VITVI_PN40024_12x_v2_chroms_URGI.tmp.fa VITVI_PN40024_12x_v2_chroms_URGI  >& bowtie2-build_VITVI_PN40024_12x_v2_chroms_URGI.log
rm VITVI_PN40024_12x_v2_chroms_URGI.tmp.fa

echo "done!"
date
