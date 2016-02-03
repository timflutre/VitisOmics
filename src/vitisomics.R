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

library(Biostrings)

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
## task: convert SNP data of the 18K Illumina array from xls to csv.gz

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

## ---------------------------------------------------------------------------
## task: extract Illumina SNP array probes into fasta files

library(Biostrings)

f <- "results/urgi/GrapeReSeq_SNP_table_180413.txt.gz"
snp.table <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)

f <- "results/urgi/GrapeReSeq_Illumina_20K_SNP_chip.txt.gz"
dat <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
table(dat$Source)
## ICVV IGA-MAF0.05  IGA-MAF0.1  URGI
##  205        1470       13347  4978

## "dat" contains too many sequences
nrow(snp.table) # 18071
nrow(dat) # 20000

sum(snp.table$Name %in% dat$Locus_Name) # 18071
dat <- dat[dat$Locus_Name %in% snp.table$Name,]
dim(dat) # 18071
table(dat$Source)
## ICVV IGA-MAF0.05  IGA-MAF0.1  URGI
##  187        1335       12040  4509

## make the sequences keeping the first allele
tmp <- head(dat)
tmp$Sequence[1]
strsplit(tmp$Sequence[1], "\\[|\\/|\\]")[[1]]
paste(strsplit(tmp$Sequence[1], "\\[|\\/|\\]")[[1]][c(1,2,4)], collapse="")
dat$Sequence.valid <- sapply(strsplit(dat$Sequence, "\\[|\\/|\\]"),
                             function(x){
                               paste(x[c(1,2,4)], collapse="")
                             })

## save valid sequences in fasta files
seq.set <- DNAStringSet(x=setNames(dat$Sequence.valid,
                                   dat$Locus_Name))
out.file <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_probes.fa.gz"
writeXStringSet(seq.set, out.file, compress=TRUE, format="fasta", width=60)

out.file <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_vinifera_probes.fa.gz"
writeXStringSet(seq.set[dat$Source != "URGI"],
                out.file, compress=TRUE, format="fasta", width=60)

out.file <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_species_probes.fa.gz"
writeXStringSet(seq.set[dat$Source == "URGI"],
                out.file, compress=TRUE, format="fasta", width=60)

## extract both alleles and their location within the sequences
dat$allele1 <- sapply(strsplit(dat$Sequence, "\\[|\\/|\\]"), `[`, 2)
dat$allele2 <- sapply(strsplit(dat$Sequence, "\\[|\\/|\\]"), `[`, 3)
dat$snp.coord.intra <- sapply(strsplit(dat$Sequence, "\\[.*"), nchar) + 1

## save information about the 18K SNPs in a file
out.file <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_array.txt.gz"
write.table(dat, out.file, quote=FALSE, sep="\t", row.names=FALSE)

## ---------------------------------------------------------------------------
## task: check the coordinates of the Illumina SNP array probes

## load alignments of "vinifera" probes on 12x0 with megablast
f <- "results/urgi/Ill18Kprobes-vinifera_12x0-chroms_megablast.txt.gz"
alns <- read.table(f, col.names=c("qseqid", "sseqid", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend",
                                  "sstart", "send", "evalue", "bitscore"),
                   stringsAsFactors=FALSE)
nrow(alns) # 15357

## some probes have multiple alignments (keep the best)
anyDuplicated(alns$qseqid) # 40
alns <- alns[! duplicated(alns$qseqid),]
summary(alns$pident) # min=89 q1=99 med=99 mean=99 q3=100 max=100
summary(alns$length) # min=62 q1=101 med=101 mean=102 q3=101 max=628

## load info about "vinifera" probes
f <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_array.txt.gz"
dat <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat <- dat[dat$Source != "URGI",]
nrow(dat) # 13562
summary(sapply(dat$Sequence.valid, nchar)) # min=84 q1=101 med=101 mean=102 q3=101 max=628

## some probes don't have any alignment
sum(! dat$Locus_Name %in% alns$qseqid) # 33

## check alignment coordinates
alns <- alns[order(alns$qseqid),]
dat2 <- dat[dat$Locus_Name %in% alns$qseqid,]
dat2 <- dat2[order(dat2$Locus_Name),]
dat2$aln.sseqid <- alns$sseqid
dat2$aln.sstart <- alns$sstart
dat2$aln.send <- alns$send
dat2$aln.pident <- alns$pident
dat2$aln.len <- alns$length

## some probes are aligned on different chromosomes than indicated
## most of them are said to belong to plastid genomes
## they were all designed by the ICVV
sum(dat2$aln.sseqid != dat2$Chromosome) # 24
dat2[dat2$aln.sseqid != dat2$Chromosome, -grep("Sequence",colnames(dat2))]

## among the probes aligned on the indicated chromosome, the indicated SNP
## coordinate is inside the alignment boundaries for all of them
sum(dat2[dat2$aln.sseqid == dat2$Chromosome, "aln.sstart"] >=
    dat2[dat2$aln.sseqid == dat2$Chromosome, "Coordinate"] ||
    dat2[dat2$aln.sseqid == dat2$Chromosome, "aln.send"] <=
    dat2[dat2$aln.sseqid == dat2$Chromosome, "Coordinate"]) # 0
