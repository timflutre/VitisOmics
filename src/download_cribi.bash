#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# http://genomes.cribi.unipd.it/DATA/

wget --timestamping http://genomes.cribi.unipd.it/DATA/README.txt

wget --timestamping http://genomes.cribi.unipd.it/DATA/Gene_V1_file_conversion.txt
wget --timestamping http://genomes.cribi.unipd.it/DATA/gene_conversion_V1_V0.txt
wget --timestamping http://genomes.cribi.unipd.it/DATA/Matching_between_V1_V2.xls

mkdir -p V2
cd V2/

wget --timestamping http://genomes.cribi.unipd.it/DATA/V2/README
wget --timestamping http://genomes.cribi.unipd.it/DATA/V2/changeLog.txt

mkdir -p V2.1
cd V2.1/

if [ ! -f V2.1.gff3.gz ]; then
  wget --timestamping http://genomes.cribi.unipd.it/DATA/V2/V2.1/V2.1.gff3
  gzip V2.1.gff3
fi

cd .. # up of V2.1/

cd .. # up of V2/
