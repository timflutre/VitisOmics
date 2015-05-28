#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# http://www.htslib.org/doc/#manual-pages
# http://broadinstitute.github.io/picard/index.html

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

date
echo "start ..."

if [ ! -f VITVI_PN40024_12x_v2_chroms_URGI.fa ]; then
  if [ ! -f VITVI_PN40024_12x_v2_chroms_URGI.fa.gz ]; then
    echo -e "ERROR: missing file VITVI_PN40024_12x_v2_chroms_URGI.fa.gz\n" 1>&2
    exit 1
  fi
  zcat VITVI_PN40024_12x_v2_chroms_URGI.fa.gz > VITVI_PN40024_12x_v2_chroms_URGI.fa
fi

samtools faidx VITVI_PN40024_12x_v2_chroms_URGI.fa >& samtools_faidx_VITVI_PN40024_12x_v2_chroms_URGI.log

java -Xmx4g -jar `which picard.jar` CreateSequenceDictionary REFERENCE=VITVI_PN40024_12x_v2_chroms_URGI.fa.gz OUTPUT=VITVI_PN40024_12x_v2_chroms_URGI.dict GENOME_ASSEMBLY=12x_v2 SPECIES="Vitis vinifera" TRUNCATE_NAMES_AT_WHITESPACE=true >& picard-CreateSequenceDictionary_VITVI_PN40024_12x_v2_chroms_URGI.log

echo "done!"
date
