## Aim: compare original files from URGI and NCBI
## Copyright (C) 2015-2016 Institut National de la Recherche Agronomique
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

## ---------------------------------------------------------------------------
## task: convert SNP data of the 18K Illumina chip from xls to csv.gz

options(java.parameters="-Xmx1024m")
library(XLConnect)

in.file <- "data/urgi/GrapeReSeq_SNP Table_all_180413.xlsx"

snp.table <- readWorksheetFromFile(file=in.file, sheet="SNP_TABLE")
## takes ~ 10 seconds

str(snp.table) # 18071 rows and 41 columns
unique(snp.table$Chr) # 34
unique(snp.table$Aux) # explained in the other sheet named "READ ME"

readme <- readWorksheetFromFile(file=in.file, sheet="READ ME")

out.file <- "results/urgi/GrapeReSeq_SNP_table_180413.txt.gz"

## write the content of the "READ ME" sheet as comments
txt <- paste0("# ", gsub("\\.", " ", colnames(readme)[1]))
cat(txt, file=gzfile(out.file), append=FALSE)
for(i in 1:nrow(readme)){
  txt <- paste0("# ", readme[i,])
  cat(txt, file=gzfile(out.file), append=TRUE)
}

## append the content of the "SNP_TABLE" sheet
write.table(x=snp.table, file=gzfile(out.file), append=TRUE, quote=FALSE,
            sep="\t", row.names=FALSE, col.names=TRUE)
## ignore the warning

in.file <- "data/urgi/GrapeReSeq_Illumina_20K_SNP_chip.xls"

dat <- readWorksheetFromFile(file=in.file, sheet="Feuil1", startRow=18)

str(dat) # 20000 rows and 20 columns
unique(dat$Genome_Build_Version) # 12 and 12.1 ?!
unique(dat$Chromosome) # 34
unique(dat$Source) # IGA-MAF0.1 IGA-MAF0.05 ICVV URGI
length(unique(dat$Ilmn_Id)) # 20000

out.file <- "results/urgi/GrapeReSeq_Illumina_20K_SNP_chip.txt.gz"
write.table(x=dat, file=gzfile(out.file), quote=FALSE, sep="\t",
            row.names=FALSE, col.names=TRUE)
