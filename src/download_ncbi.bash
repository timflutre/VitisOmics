#!/usr/bin/env bash

# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/README
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/README_CURRENT_BUILD

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/allcontig.agp.gz

mkdir -p Assembled_chromosomes/; cd Assembled_chromosomes/
for i in {1..19}; do
    wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/Assembled_chromosomes/seq/vvi_ref_12X_chr${i}.fa.gz
done
for x in {"chrMT","chrPltd","unlocalized","unplaced"}; do
    wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/Assembled_chromosomes/seq/vvi_ref_12X_${x}.fa.gz
done
cd ..

mkdir -p GFF; cd GFF
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/GFF/ref_12X_scaffolds.gff3.gz
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/GFF/ref_12X_top_level.gff3.gz
cd ..
