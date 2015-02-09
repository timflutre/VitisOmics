#!/usr/bin/env bash

# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/doc/blast/formatdb.html

if [ ! -f Vvin-PN40024-12x-chr.fa.gz ]; then
    echo -e "ERROR: can't find file Vvin-PN40024-12x-chr.fa.gz\n" 1>&2
    exit 1
fi

zcat Vvin-PN40024-12x-chr.fa.gz \
    | formatdb -i stdin -l formatdb_Vvin-PN40024-12x-chr.log -p F -n Vvin-PN40024-12x-chr
