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
## task: compare URGI and NCBI fasta files for 12x scaffolds

library(Biostrings) # http://bioconductor.org

urgi <-
  readDNAStringSet(filepath="results/urgi/VITVI_PN40024_12x_v0_scaffolds_EMBL_r102.fa.gz",
                   format="fasta")
length(urgi) # 2059
names(urgi[1]) # "scaffold_0 FN594950|Vitis vinifera|PN40024|assembly12x.0"
width(urgi[1]) # 13101952
letterFrequency(urgi[1], letters=c("A","T","G","C","N"))

ncbi <-
  readDNAStringSet(filepath="results/ncbi/VITVI_PN40024_12x_v0_scaffolds_NCBI.fa.gz",
                   format="fasta")
length(ncbi) # 2059
names(ncbi[1]) # "scaffold_11 NW_003723981.1|Vitis vinifera|PN40024|assembly12x.0|localized chr1"
width(ncbi[1]) # 6551432
letterFrequency(ncbi[1], letters=c("A","T","G","C","N"))

fa.head <- data.frame(urgi=sapply(strsplit(names(urgi), " "), function(x){x[1]}),
                      ncbi=sapply(strsplit(names(ncbi), " "), function(x){x[1]}),
                      stringsAsFactors=FALSE)

all(sort(fa.head$urgi) == sort(fa.head$ncbi)) # TRUE

for(i in 1:nrow(fa.head)){
  if(urgi[which(fa.head$urgi == fa.head$urgi[i])] !=
     ncbi[which(fa.head$ncbi == fa.head$urgi[i])]){
    message(fa.head$urgi[i])
    break
  }
} ## all good!
