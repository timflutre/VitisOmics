#!/usr/bin/env bash

# Copyright (C) 2014-2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/README
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/README_CURRENT_BUILD

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/scaffold_names
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/allcontig.agp.gz

mkdir -p Assembled_chromosomes; cd Assembled_chromosomes/
for i in {1..19}; do
  wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/Assembled_chromosomes/seq/vvi_ref_12X_chr${i}.fa.gz
done
for x in {"chrMT","chrPltd","unlocalized","unplaced"}; do
  wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/Assembled_chromosomes/seq/vvi_ref_12X_${x}.fa.gz
done
cd ..

mkdir -p GFF; cd GFF/
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/GFF/ref_12X_scaffolds.gff3.gz
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/GFF/ref_12X_top_level.gff3.gz
cd ..

mkdir -p ARCHIVE
cd ARCHIVE/
mkdir -p BUILD.1.1
cd BUILD.1.1/

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/README
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/README_CURRENT_BUILD

wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/scaffold_names

mkdir -p Assembled_chromosomes
cd Assembled_chromosomes/
for i in {1..19}; do
  wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/Assembled_chromosomes/vvi_ref_chr${i}.fa.gz
done
cd ..

mkdir -p CHRS; cd CHRS/
for i in {1..9}; do
  wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/CHR_0${i}/vvi_ref_chr${i}.fa.gz
done
for i in {10..19}; do
  wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/CHR_${i}/vvi_ref_chr${i}.fa.gz
done
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/CHR_Pltd/vvi_ref_chrPltd.fa.gz
wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/Vitis_vinifera/ARCHIVE/BUILD.1.1/CHR_Un/vvi_ref_chrUn.fa.gz
cd ../

cd ../..
