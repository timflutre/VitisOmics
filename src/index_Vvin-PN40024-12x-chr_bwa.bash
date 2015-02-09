#!/usr/bin/env bash

# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# http://bio-bwa.sourceforge.net/bwa.shtml

if [ ! -f Vvin-PN40024-12x-chr.fa.gz ]; then
  echo -e "ERROR: missing file Vvin-PN40024-12x-chr.fa.gz\n" 1>&2
  exit 1
fi

bwa index -p Vvin-PN40024-12x-chr Vvin-PN40024-12x-chr.fa.gz >& bwa_index_Vvin-PN40024-12x-chr.log
