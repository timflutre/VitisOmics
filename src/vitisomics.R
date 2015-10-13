## Aim: compare original files from URGI and NCBI
## Copyright (C) 2015 Institut National de la Recherche Agronomique
## License: GPL-3+
## Persons: Timothée Flutre [cre,aut]
## Versioning: https://github.com/timflutre/VitisOmics/src

rm(list=ls())

## path to the VitisOmics repository on your computer
repo.dir <- "<...>/VitisOmics"
setwd(repo.dir)

## ---------------------------------------------------------------------------
## task: look at the AGP file for "PN40024 8x" from NCBI

## read if necessary: http://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/

## load sequences
cmd <- paste0("zcat ", repo.dir,
              "/data/ncbi/ARCHIVE/BUILD.1.1/allcontig.agp.gz",
              " | sed -n 1p | sed 's/#//g'",
                "; zcat ", repo.dir,
              "/data/ncbi/ARCHIVE/BUILD.1.1/allcontig.agp.gz",
              " | grep -v \"#\" | awk '{if($5==\"W\")print $0}'")
agp.seq <- read.table(file=pipe(cmd), header=TRUE)
str(agp.seq) # 19577 rows

## load gaps
cmd <- paste0("zcat ", repo.dir,
              "/data/ncbi/ARCHIVE/BUILD.1.1/allcontig.agp.gz",
              " | sed -n 2p",
              " | awk '{printf $3; for(i=4;i<=NF;++i)printf \"\\t\"$i; printf \"\\n\" }'",
              "; zcat ", repo.dir,
              "/data/ncbi/ARCHIVE/BUILD.1.1/allcontig.agp.gz",
              " | grep -v \"#\" | awk '{if($5==\"N\")print $0}'")
agp.gap <- read.table(file=pipe(cmd), header=TRUE, sep="\t")
str(agp.gap) # 16063 rows

## ---------------------------------------------------------------------------
## task: compare URGI and NCBI fasta file

library(Biostrings) # http://bioconductor.org

VITVI.PN40024.8x.chroms.URGI <-
  readDNAStringSet(filepath="results/urgi/VITVI_PN40024_8x_chroms_URGI.fa.gz",
                   format="fasta")
length(VITVI.PN40024.8x.chroms.URGI) # 35
names(VITVI.PN40024.8x.chroms.URGI) # chr1, chr10, chr10_random, etc
width(VITVI.PN40024.8x.chroms.URGI)
letterFrequency(VITVI.PN40024.8x.chroms.URGI, letters=c("A","T","G","C","N"))

VITVI.PN40024.8x.chroms.NCBI <-
  readDNAStringSet(filepath="results/ncbi/VITVI_PN40024_8x_chroms_NCBI.fa.gz",
                   format="fasta")
length(VITVI.PN40024.8x.chroms.NCBI) # 3343
head(names(VITVI.PN40024.8x.chroms.NCBI), n=20) #
width(VITVI.PN40024.8x.chroms.NCBI)
letterFrequency(VITVI.PN40024.8x.chroms.NCBI, letters=c("A","T","G","C","N"))

names(VITVI.PN40024.8x.chroms.URGI[1])
names(VITVI.PN40024.8x.chroms.NCBI[1])
VITVI.PN40024.8x.chroms.URGI[1] == VITVI.PN40024.8x.chroms.NCBI[1] # TRUE
