#!/usr/bin/env bash

# Copyright (C) 2015-2016 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# https://urgi.versailles.inra.fr/Species/Vitis/Data-Sequences/Genome-sequences
# https://urgi.versailles.inra.fr/Species/Vitis/Annotations

echo "12X.0 version of the grapevine reference genome sequence from The French-Italian Public Consortium (PN40024)"
if [ ! -f VV_12X_embl_102_WGS_contigs.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_12X_embl_102_WGS_contigs.fsa.zip
fi
if [ ! -f VV_12X_embl_102_Scaffolds.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_12X_embl_102_Scaffolds.fsa.zip
fi
if [ ! -f VV_chr12x.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_chr12x.fsa.zip
fi
if [ $(ls -1 *.zip 2>/dev/null | wc -l) -ne 0 ]; then
  ls *.zip | while read f; do
    if [ ! -f "${prefix}.gz" ]; then
      unzip $f
      prefix=$(echo $f | sed 's/.zip//')
      gzip $prefix
      rm -f $f
    fi
  done
fi

echo "12X.0 \"golden path\" files for chromosome assembly (2009_12_04)"
## note: wget doesn't allow "--timestamping" and "-O" together
if [ ! -f 12x0_chr.agp ]; then
  wget --timestamping https://urgi.versailles.inra.fr/content/download/1028/8244/file/chr.agp
  mv chr.agp 12x0_chr.agp
fi
if [ ! -f 12x0_chrUn.agp ]; then
  wget --timestamping https://urgi.versailles.inra.fr/content/download/1029/8248/file/chrUn.agp
  mv chrUn.agp 12x0_chrUn.agp
fi
if [ ! -f 12x0_chr.agp.info ]; then
  wget --timestamping https://urgi.versailles.inra.fr/content/download/2149/19329/file/chr.agp.info
  mv chr.agp.info 12x0_chr.agp.info
fi
if [ ! -f 12x0_chr.lg ]; then
  wget --timestamping https://urgi.versailles.inra.fr/content/download/2150/19333/file/chr.lg
  mv chr.lg 12x0_chr.lg
fi
if [ ! -f 12x0_scaffolds.lg ]; then
  wget --timestamping https://urgi.versailles.inra.fr/content/download/1093/8684/file/scaffolds.lg
  mv scaffolds.lg 12x0_scaffolds.lg
fi

echo "12X.2 version of the grapevine reference genome sequence from The French-Italian Public Consortium (PN40024)"
wget --timestamping https://urgi.versailles.inra.fr/download/vitis/12Xv2_grapevine_genome_assembly.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/3044/26115/file/golden_path_V2_111113_allChr.csv
wget --timestamping https://urgi.versailles.inra.fr/content/download/4702/36526/file/vitis12XV2_final.agp
wget --timestamping https://urgi.versailles.inra.fr/content/download/3043/26111/file/chr_size_V2.txt

echo "V0 annotation of the 12X.0 genome assembly (Genoscope)"
mkdir -p 12x_annotation_Genoscope_V0; cd 12x_annotation_Genoscope_V0
wget --timestamping https://urgi.versailles.inra.fr/content/download/2157/19376/file/Vitis_vinifera_annotation.gff.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2158/19380/file/Vitis_vinifera_mRNA.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2159/19384/file/Vitis_vinifera_peptide.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2199/19692/file/RepeatMasker.gff3.bz2
wget --timestamping https://urgi.versailles.inra.fr/content/download/2200/19696/file/TRF.gff3.bz2
cd ..

echo "V1 annotation of the 12X.0 genome assembly (CRIBI)"
mkdir -p 12x_annotation_CRIBI_V1; cd 12x_annotation_CRIBI_V1/
wget --timestamping https://urgi.versailles.inra.fr/content/download/2160/19388/file/chrAll.jigsawgaze_NR.gff.gz
https://urgi.versailles.inra.fr/content/download/2161/19392/file/chrAll.jigsawgaze_RO.gff.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2163/19400/file/V1_JigsawGazeNR_mRNA.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2164/19404/file/V1_JigsawGazeRO_mRNA.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2166/19412/file/V1_JigsawGazeNR.prot.fa.gz
wget --timestamping https://urgi.versailles.inra.fr/content/download/2167/19416/file/V1_JigsawGazeRO.prot.fa.gz
cd ..

echo "8X version of the of the grapevine reference genome sequence from The French-Italian Public Consortium (PN40024)"
if [ ! -f VV_chr8x.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_chr8x.fsa.zip
  unzip VV_chr8x.fsa.zip
  gzip VV_chr8x.fsa
  rm -f VV_chr8x.fsa.zip
fi
if [ ! -f VV_8X_embl_98_Scaffolds.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_8X_embl_98_Scaffolds.fsa.zip
  unzip VV_8X_embl_98_Scaffolds.fsa.zip
  gzip VV_8X_embl_98_Scaffolds.fsa
  rm -f VV_8X_embl_98_Scaffolds.fsa.zip
fi
if [ ! -f VV_8X_embl_98_WGS_contigs.fsa.gz ]; then
  wget --timestamping https://urgi.versailles.inra.fr/download/vitis/VV_8X_embl_98_WGS_contigs.fsa.zip
  unzip VV_8X_embl_98_WGS_contigs.fsa.zip
  gzip VV_8X_embl_98_WGS_contigs.fsa
  rm -f VV_8X_embl_98_WGS_contigs.fsa.zip
fi

echo "Illumina Infinium chip (GrapeReSeq, 18K)"
wget --timestamping https://urgi.versailles.inra.fr/content/download/2688/23435/file/GrapeReSeq_SNP%20Table_all_180413.xlsx
wget --timestamping https://urgi.versailles.inra.fr/content/download/2689/23439/file/GrapeReSeq_CF-EPGV-v2b.egt
wget --timestamping https://urgi.versailles.inra.fr/content/download/2465/21587/file/GrapeReSeq_Illumina_20K_SNP_chip.xls
