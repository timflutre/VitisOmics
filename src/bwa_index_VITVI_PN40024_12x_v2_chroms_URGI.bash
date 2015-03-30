#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]

# http://bio-bwa.sourceforge.net/bwa.shtml

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

date
echo "start ..."

if [ ! -f VITVI_PN40024_12x_v2_chroms_URGI.fa.gz ]; then
  echo -e "ERROR: missing file VITVI_PN40024_12x_v2_chroms_URGI.fa.gz\n" 1>&2
  exit 1
fi

bwa index -p VITVI_PN40024_12x_v2_chroms_URGI VITVI_PN40024_12x_v2_chroms_URGI.fa.gz >& bwa_index_VITVI_PN40024_12x_v2_chroms_URGI.log

echo "done!"
date
