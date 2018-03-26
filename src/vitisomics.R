## Aim: compare original files from URGI and NCBI
## Copyright (C) 2015-2018 Institut National de la Recherche Agronomique
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
## task: compare URGI and CRIBI fasta files for 12x.0 chromosomes

library(Biostrings)

## load URGI file
urgi <-
  readDNAStringSet(filepath="data/urgi/VV_chr12x.fsa.gz",
                   format="fasta")
length(urgi) # 33
names(urgi[1]) # "FN597015 Vitis vinifera, line PN40024, unoriented chromosome_1, chr1 "
width(urgi[1]) # 23037639
letterFrequency(urgi[1], letters=c("A","T","G","C","N"))

## load CRIBI file
tmp.file <- "data/cribi/Genome12X_all-chrs.fa.gz"
if(! file.exists(tmp.file)){
  cmd <- paste0("tar -xzOf data/cribi/Genome12X.tar.gz | gzip > ", tmp.file)
  system(cmd)
}
cribi <- readDNAStringSet(filepath=tmp.file, format="fasta")
length(cribi) # 33
names(cribi[1]) # "chr1"
width(cribi[1]) # 23037639
letterFrequency(cribi[1], letters=c("A","T","G","C","N"))

fa.head <- data.frame(urgi=sapply(strsplit(names(urgi), " "), function(x){x[length(x)]}),
                      cribi=names(cribi),
                      stringsAsFactors=FALSE)

all(sort(fa.head$urgi) == sort(fa.head$cribi)) # TRUE

for(i in 1:nrow(fa.head)){
  if(urgi[which(fa.head$urgi == fa.head$urgi[i])] !=
     cribi[which(fa.head$cribi == fa.head$urgi[i])]){
    message(fa.head$urgi[i])
    break
  }
} ## all good!

## ---------------------------------------------------------------------------
## task: convert SNP data of the 18K Illumina array from xls to txt.gz

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
write.table(dat, gzfile(out.file), quote=FALSE, sep="\t", row.names=FALSE)

## ---------------------------------------------------------------------------
## task: plot SNP density along 12x0 chromosomes from Illumina 18k array

library(ggbio)

## load info about the SNPs
f <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_array.txt.gz"
dat <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
str(dat) # 18071 x 24

## keep only "vinifera" SNPs on non-random chromosomes
dat <- dat[dat$Source != "URGI",]
nrow(dat) # 13562
table(dat$Chromosome)
dat <- dat[! grepl("random", dat$Chromosome),]
dat <- dat[! grepl("chrUn", dat$Chromosome),]
dat <- dat[! grepl("Pltd", dat$Chromosome),]
nrow(dat) # 11853

## get "seqinfo" from 12x0
library(BSgenome.Vvinifera.URGI.IGGP12Xv0)
seqinfo12x0 <- seqinfo(BSgenome.Vvinifera.URGI.IGGP12Xv0)

## convert SNP coords into GRanges
gr <- GRanges(seqnames=Rle(dat$Chromosome),
              ranges=IRanges(start=dat$Coordinate, end=dat$Coordinate),
              strand=rep("*", nrow(dat)),
              seqinfo=seqinfo12x0)
seqlevels(gr) <- seqlevelsInUse(gr)
gr
seqlevels(seqinfo12x0) <- seqlevels(gr)

## plot the karyogram with the SNP density
autoplot(seqinfo12x0, layout="karyogram", main="Karyogram of Vitis vinifera L.")

## plot the SNP density over the karyogram
autoplot(gr, layout="karyogram",
         main=paste0("SNP density (", length(gr), " SNPs"))

## ---------------------------------------------------------------------------
## task: check the alignment coordinates of the Illumina SNP array probes

## load alignments of "vinifera" probes with megablast
min.ver <- 2
f <- paste0("results/urgi/Ill18Kprobes-vinifera_12x", min.ver,
            "-chroms_megablast.txt.gz")
alns <- read.table(f, col.names=c("qseqid", "sseqid", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend",
                                  "sstart", "send", "evalue", "bitscore"),
                   stringsAsFactors=FALSE)
nrow(alns) # min.ver=0: 15357 ; min.ver=2: 15358

## some probes have multiple alignments (keep the best)
anyDuplicated(alns$qseqid) # 40 for both
alns <- alns[! duplicated(alns$qseqid),]
summary(alns$pident) # min.ver=0: min=89 q1=99 med=99 mean=99 q3=100 max=100
summary(alns$length) # min.ver=0: min=62 q1=101 med=101 mean=102 q3=101 max=628

## load info about "vinifera" probes
f <- "results/urgi/GrapeReSeq_Illumina_18K_SNP_array.txt.gz"
dat <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat <- dat[dat$Source != "URGI",]
nrow(dat) # 13562
summary(sapply(dat$Sequence.valid, nchar)) # min=84 q1=101 med=101 mean=102 q3=101 max=628

## some probes don't have any alignment
sum(! dat$Locus_Name %in% alns$qseqid) # 33 for both

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
sum(dat2$aln.sseqid != dat2$Chromosome) # min.ver=0: 24 ; min.ver=2: 1723 ("random" chr...)
head(dat2[dat2$aln.sseqid != dat2$Chromosome, -grep("Sequence",colnames(dat2))])

## among the probes aligned on the indicated chromosome, the indicated SNP
## coordinate is inside the alignment boundaries for all of them
sum(dat2[dat2$aln.sseqid == dat2$Chromosome, "aln.sstart"] >=
    dat2[dat2$aln.sseqid == dat2$Chromosome, "Coordinate"] ||
    dat2[dat2$aln.sseqid == dat2$Chromosome, "aln.send"] <=
    dat2[dat2$aln.sseqid == dat2$Chromosome, "Coordinate"])
## min.ver=0: 0 ; min.ver=2: 1

## ---------------------------------------------------------------------------
## task: extract Illumina SNP metadata on 12x0 and 12x2

options(java.parameters="-Xmx1024m")
library(XLConnect)

setwd(paste0(repo.dir, "/results/grapereseq_18k_vitis_microarray"))

## extract the xlsx file
p2f <- paste0(repo.dir, "/data/urgi/Grapereseq_cultivated_pool_data.zip")
(files <- unzip(zipfile=p2f, list=TRUE))
p2f.xlsx <- files$Name[grep("xlsx", files$Name)]
if(file.exists(basename(p2f.xlsx)))
  file.remove(basename(p2f.xlsx))
unzip(zipfile=p2f, files=p2f.xlsx, junkpaths=TRUE)

## convert xlsx into tsv
system.time(
    dat <- readWorksheetFromFile(file=basename(p2f.xlsx),
                                 sheet="SNP Markers"))
## takes ~ 6 seconds
str(dat)
colnames(dat)[colnames(dat) == "Variation.Allele_réel"] <- "Variation.RealAllele"
colnames(dat)[colnames(dat) == "X5.Flanking"] <- "Flanking.5"
colnames(dat)[colnames(dat) == "X3.Flanking"] <- "Flanking.3"

## write in gzip-compressed tsv format
tsv.file <- "grapereseq_18k_vitis_microarray_12x0-2-12x2.tsv.gz"
write.table(dat, file=gzfile(tsv.file), quote=FALSE, sep="\t", row.names=FALSE)

## ---------------------------------------------------------------------------
## task: make TxDb on IGGP12Xv0 from CRIBI (V2.1)

library(rtracklayer)

setwd("results/make_TxDb_IGGP12Xv0_CRIBIv2-1")

out.dir <- "input_AnnotHub_IGGP12Xv0_CRIBIv2.1"
if(dir.exists(out.dir))
  unlink(out.dir, recursive=TRUE)
dir.create(path=out.dir)

p2f <- "V2.1_updated.gff3.gz"
file.copy(from=p2f, to=paste0(out.dir, "/", basename(p2f)))

gr <- import.gff(con=p2f, version="3", genome="IGGP12Xv0",
                 sequenceRegionsAsSeqinfo=TRUE)
length(gr) # 820944

p2f <- paste0(out.dir, "/GRanges.RData")
save(gr, file=p2f)

p2f <- paste0(out.dir, "/metadata.txt")
if(file.exists(p2f))
  file.remove(p2f)
cat("SourceUrl: http://genomes.cribi.unipd.it/DATA/V2/V2.1/V2.1.gff3\n",
    file=p2f, append=TRUE)
cat("SourceType: gff3\n", file=p2f, append=TRUE)
cat("SourceVersion: 2.1\n", file=p2f, append=TRUE)
cat("SourceLastModifiedDate: 2014-04-17\n", file=p2f, append=TRUE)
cat("DataProvider: CRIBI\n", file=p2f, append=TRUE)
cat("Title: Vvinifera_CRIBI_IGGP12Xv0_V2.1.gff3.Rdata\n", file=p2f, append=TRUE)
cat("Description: Gene Annotation for Vitis vinifera\n", file=p2f, append=TRUE)
cat("Species: Vitis vinifera\n", file=p2f, append=TRUE)
cat("TaxonomyId: 29760\n", file=p2f, append=TRUE)
cat("Genome: IGGP12Xv0\n", file=p2f, append=TRUE)
cat("Maintainer: Timothée Flutre timothee.flutre@supagro.inra.fr\n", file=p2f, append=TRUE)
cat("Notes: compare to the original file, meta-data were added, and chrUn was renamed into chrUkn\n", file=p2f, append=TRUE)

tar(tarfile=paste0(out.dir, ".tar.gz"),
    files=out.dir, compression="gzip")
## to be sent to Bioconductor

## ---------------------------------------------------------------------------
## task: check the TxDb on IGGP12Xv0 from CRIBI (V2.1)

library(AnnotationHub)

ahub <- AnnotationHub()
"Vitis vinifera" %in% unique(ahub$species)
"CRIBI" %in% unique(ahub$dataprovider)
query(ahub, c("Vitis vinifera", "CRIBI"))

## download the GRanges
gr <- ahub[["AH50773"]]

## make the TxDb
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gr)

## save the TxDb into a ".sqlite" database file
## so that it can be made available to other users
p2f <- paste0("results/make_TxDb_IGGP12Xv0_CRIBIv2-1/",
              "TxDb_Vvinifera_IGGP12Xv0_CRIBIv2-1.sqlite")
saveDb(x=txdb, file=p2f)

## load the TxDb
txdb <- loadDb(file=p2f)

## have a look at the resource
txdb
length(genes(txdb)) # 31845
length(transcripts(txdb)) # 55564
length(exons(txdb)) # 321050
length(cds(txdb)) # 297312

## let us choose a "good-example" gene:
## VIT_201s0011g00050: chr1, 3 mRNAs, 13 exons
gene.name <- "VIT_201s0011g00050"

g <- genes(txdb)
g[gene.name]

t <- transcriptsBy(txdb, "gene")
t[gene.name]

eg <- exonsBy(txdb, "gene")
eg[gene.name]

et <- exonsBy(txdb, "tx", use.names=TRUE)
et[paste0(gene.name, ".2")]

f <- fiveUTRsByTranscript(txdb, use.names=TRUE)
f[paste0(gene.name, ".2")] # note that the 5' UTR IDs from the GFF3 file are absent

## ---------------------------------------------------------------------------
## task: make TxDb on IGGP12Xv0 from Genoscope

library(rtracklayer)

setwd("results/make_TxDb_IGGP12Xv0_Genoscope")

out.dir <- "input_AnnotHub_IGGP12Xv0_Genoscope"
if(dir.exists(out.dir))
  unlink(out.dir, recursive=TRUE)
dir.create(path=out.dir)

p2f <- "Vitis_vinifera_annotation_updated.gff3.gz"
file.copy(from=p2f, to=paste0(out.dir, "/", basename(p2f)))

gr <- import.gff(con=p2f, version="3", genome="IGGP12Xv0",
                 sequenceRegionsAsSeqinfo=TRUE)
length(gr) # 209457

p2f <- paste0(out.dir, "/GRanges.RData")
save(gr, file=p2f)

p2f <- paste0(out.dir, "/metadata.txt")
if(file.exists(p2f))
  file.remove(p2f)
cat("SourceUrl: http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/annotation/Vitis_vinifera_annotation.gff.gz\n",
    file=p2f, append=TRUE)
cat("SourceType: gff3\n", file=p2f, append=TRUE)
cat("SourceVersion: 1.0\n", file=p2f, append=TRUE)
cat("SourceLastModifiedDate: 2010-03-19\n", file=p2f, append=TRUE)
cat("DataProvider: Genoscope\n", file=p2f, append=TRUE)
cat("Title: Vvinifera_Genoscope_IGGP12Xv0_V1.0.gff3.Rdata\n", file=p2f, append=TRUE)
cat("Description: Gene Annotation for Vitis vinifera\n", file=p2f, append=TRUE)
cat("Species: Vitis vinifera\n", file=p2f, append=TRUE)
cat("TaxonomyId: 29760\n", file=p2f, append=TRUE)
cat("Genome: IGGP12Xv0\n", file=p2f, append=TRUE)
cat("Maintainer: Timothée Flutre timothee.flutre@supagro.inra.fr\n", file=p2f, append=TRUE)
cat("Notes: compare to the original file, the format was upgraded from GFF2 to GFF3 only keeping rows corresponding to gene/mRNA/CDS, meta-data were added, and chrUn was renamed into chrUkn\n", file=p2f, append=TRUE)

tar(tarfile=paste0(out.dir, ".tar.gz"),
    files=out.dir, compression="gzip")
## to be sent to Bioconductor

## ---------------------------------------------------------------------------
## task: check the TxDb on IGGP12Xv0 from Genoscope

library(AnnotationHub)

ahub <- AnnotationHub()
"Vitis vinifera" %in% unique(ahub$species)
"Genoscope" %in% unique(ahub$dataprovider)
query(ahub, c("Vitis vinifera", "Genoscope"))

## download the GRanges
gr <- ahub[["AH50774"]]

## make the TxDb
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gr)

## save the TxDb into a ".sqlite" database file
## so that it can be made available to other users
p2f <- paste0("results/make_TxDb_IGGP12Xv0_Genoscope/",
              "TxDb_Vvinifera_IGGP12Xv0_Genoscope.sqlite")
saveDb(x=txdb, file=p2f)

## load the TxDb
txdb <- loadDb(file=p2f)

## have a look at the resource
txdb
length(genes(txdb)) # 26346
length(transcripts(txdb)) # 26346
length(exons(txdb)) # 156765
length(cds(txdb)) # 156765

## let us choose a "good-example" gene:
## GSVIVG01000001001: chr14, 1 mRNA, 4 CDSs, 1 UTR
gene.name <- "GSVIVG01000001001"

g <- genes(txdb)
g[gene.name]

t <- transcriptsBy(txdb, "gene")
t[gene.name]

eg <- exonsBy(txdb, "gene")
eg[gene.name]

et <- exonsBy(txdb, "tx", use.names=TRUE)
et["GSVIVT01000001001"]

## ---------------------------------------------------------------------------
## task: make TxDb on IGGP8X from Genoscope

library(rtracklayer)

setwd("results/make_TxDb_IGGP8X_Genoscope")

out.dir <- "input_AnnotHub_IGGP8X_Genoscope"
if(dir.exists(out.dir))
  unlink(out.dir, recursive=TRUE)
dir.create(path=out.dir)

p2f <- "Vitis_vinifera_annotation_v1_updated.gff3.gz"
file.copy(from=p2f, to=paste0(out.dir, "/", basename(p2f)))

gr <- import.gff(con=p2f, version="3", genome="IGGP8X",
                 sequenceRegionsAsSeqinfo=TRUE)
length(gr) # 210219

p2f <- paste0(out.dir, "/GRanges.RData")
save(gr, file=p2f)

p2f <- paste0(out.dir, "/metadata.txt")
if(file.exists(p2f))
  file.remove(p2f)
cat("SourceUrl: http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/8X/annotation/Vitis_vinifera_annotation_v1.gff\n",
    file=p2f, append=TRUE)
cat("SourceType: gff3\n", file=p2f, append=TRUE)
cat("SourceVersion: 1.0\n", file=p2f, append=TRUE)
cat("SourceLastModifiedDate: 2007-10-09\n", file=p2f, append=TRUE)
cat("DataProvider: Genoscope\n", file=p2f, append=TRUE)
cat("Title: Vvinifera_Genoscope_IGGP8X_V1.0.gff3.Rdata\n", file=p2f, append=TRUE)
cat("Description: Gene Annotation for Vitis vinifera\n", file=p2f, append=TRUE)
cat("Species: Vitis vinifera\n", file=p2f, append=TRUE)
cat("TaxonomyId: 29760\n", file=p2f, append=TRUE)
cat("Genome: IGGP8X\n", file=p2f, append=TRUE)
cat("Maintainer: Timothée Flutre timothee.flutre@supagro.inra.fr\n", file=p2f, append=TRUE)
cat("Notes: compare to the original file, the format was upgraded from GFF2 to GFF3 only keeping rows corresponding to gene/mRNA/CDS, meta-data were added, and chrUn was renamed into chrUkn\n", file=p2f, append=TRUE)

tar(tarfile=paste0(out.dir, ".tar.gz"),
    files=out.dir, compression="gzip")
## to be sent to Bioconductor

## ---------------------------------------------------------------------------
## task: check the TxDb on IGGP8X from Genoscope

library(AnnotationHub)

ahub <- AnnotationHub()
"Vitis vinifera" %in% unique(ahub$species)
"Genoscope" %in% unique(ahub$dataprovider)
query(ahub, c("Vitis vinifera", "Genoscope"))

## download the GRanges
gr <- ahub[["AH50775"]]

## make the TxDb
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gr)

## save the TxDb into a ".sqlite" database file
## so that it can be made available to other users
p2f <- paste0("results/make_TxDb_IGGP8X_Genoscope/",
              "TxDb_Vvinifera_IGGP8X_Genoscope.sqlite")
saveDb(x=txdb, file=p2f)

## load the TxDb
txdb <- loadDb(file=p2f)

## have a look at the resource
txdb
length(genes(txdb)) # 30434
length(transcripts(txdb)) # 30434
length(exons(txdb)) # 149351
length(cds(txdb)) # 149351

## ---------------------------------------------------------------------------
## task: make TxDb on IGGP12Xv2 from Canaguier et al (2017) known as VCost.v3

library(rtracklayer)

setwd("results/make_TxDb_IGGP12Xv2_Canaguier2017")

out.dir <- "input_AnnotHub_IGGP12Xv2_Canaguier2017"
if(dir.exists(out.dir))
  unlink(out.dir, recursive=TRUE)
dir.create(path=out.dir)

p2f <- "VCost.v3_20.gff3.gz"
file.copy(from=p2f, to=paste0(out.dir, "/", basename(p2f)))

gr <- import.gff(con=p2f, version="3", genome="IGGP12Xv2",
                 sequenceRegionsAsSeqinfo=TRUE)
length(gr)
## with VCost.v3 (v10 from November 2017): 532145
## with VCost.v3 (v11 from December 2017): 532145
## with VCost.v3 (v14prep from December 21, 2017): 531877
## with VCost.v3 (v15 from January 22, 2018): 531877
## with VCost.v3 (v17 from January 23, 2018): 531842
## with VCost.v3 (v18 from January 24, 2018): 531841
## with VCost.v3 (v20 from February 16, 2018): 531745
sum(is.na(gr$Name)) # v17: 2; v18: 0; v20: 0

p2f <- paste0(out.dir, "/GRanges.RData")
save(gr, file=p2f)

p2f <- paste0(out.dir, "/metadata.txt")
if(file.exists(p2f))
  file.remove(p2f)
cat("SourceUrl: http://doi.org/10.15454/1.5009072354498936E12\n",
    file=p2f, append=TRUE)
cat("SourceType: gff3\n", file=p2f, append=TRUE)
cat("SourceVersion: 3.20\n", file=p2f, append=TRUE)
cat("SourceLastModifiedDate: 2018-02-16\n", file=p2f, append=TRUE)
cat("DataProvider: URGI\n", file=p2f, append=TRUE)
cat("Title: Vvinifera_URGI_IGGP12Xv2_V3-20.gff3.Rdata\n", file=p2f, append=TRUE)
cat("Description: Gene Annotation for Vitis vinifera\n", file=p2f, append=TRUE)
cat("Species: Vitis vinifera\n", file=p2f, append=TRUE)
cat("TaxonomyId: 29760\n", file=p2f, append=TRUE)
cat("Genome: IGGP12Xv2\n", file=p2f, append=TRUE)
cat("Maintainer: Timothée Flutre timothee.flutre@inra.fr\n", file=p2f, append=TRUE)
cat("Notes: compare to the original GFF3 file, chromosomes were slightly renamed to be compatible with the reference genome\n",
    file=p2f, append=TRUE)

tar(tarfile=paste0(out.dir, ".tar.gz"),
    files=out.dir, compression="gzip")
## to be sent to Bioconductor (only if the next commands have no error!)

library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gr)
## with VCost.v3 (v10 from November 2017)
##   error: 655 CDS have a missing phase
## with VCost.v3 (v11 from December 2017)
##   error: some exons are linked to transcripts not found in the file
##   warning: the following orphan exon were dropped
##   warning: the following orphan CDS were dropped
##   ...
## with VCost.v3 (v14prep from December 21, 2017)
##   no error
##   18 warnings: but only 2 seem to be related to the input data
##     dropped transcripts because their exon ranks could not be inferred
##     rejected transcripts because they have CDSs that cannot be mapped to an exon
## with VCost.v3 (v15 from January 22, 2018)
##   no error
##   2 warnings
##     dropped transcripts because their exon ranks could not be inferred
##     rejected transcripts because they have CDSs that cannot be mapped to an exon
## with VCost.v3 (v17 from January 23, 2018)
##   no error
##   2 warnings
##     dropped transcripts because their exon ranks could not be inferred
##     rejected transcripts because they have CDSs that cannot be mapped to an exon
## with VCost.v3 (v18 from January 24, 2018)
##   no error
##   2 warnings
##     dropped transcripts because their exon ranks could not be inferred
##     rejected transcripts because they have CDSs that cannot be mapped to an exon
## with VCost.v3 (v20 from February 16, 2018)
##   no error
##   1 warning
##     "phase" column contains non-NA for exon; info ignored

## ---------------------------------------------------------------------------
## task: check the TxDb on IGGP12Xv2 from Canaguier et al (2017) known as VCost.v3

library(AnnotationHub)

ahub <- AnnotationHub()
"Vitis vinifera" %in% unique(ahub$species)

query(ahub, "Vitis vinifera")

"URGI" %in% unique(ahub$dataprovider)
query(ahub, c("Vitis vinifera", "URGI"))

## download the GRanges
gr <- ahub[["AH60919"]]

## make the TxDb
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gr)

## save the TxDb into a ".sqlite" database file
## so that it can be made available to other users
p2f <- paste0("results/make_TxDb_IGGP12Xv2_Canaguier2017/",
              "TxDb_Vvinifera_IGGP12Xv2_URGIv3-20.sqlite")
if(file.exists(p2f))
  file.remove(p2f)
saveDb(x=txdb, file=p2f)

## load the TxDb
txdb <- loadDb(file=p2f)

## have a look at the resource
txdb
length(genes(txdb)) # v14prep: 42354 ; v15: 42357 ; v17: 42409; v18: 42413; v20: 42413
length(transcripts(txdb)) # v14prep: 49359 ; v15: 49359 ; v17: 49484; v18: 49487; v20: 50814
length(exons(txdb)) # v14prep: 200958 ; v15: 200958 ; v17: 201557; v18: 201595; v20: 202982
length(cds(txdb)) # v14prep: 187169 ; v15: 187169 ; v17: 227038; v18: 227124; v20: 234788

## let us choose a "good-example" gene:
## Vitvi01g00050: chr1, 3 mRNAs, 13 exons
gene.name <- "Vitvi01g00050"

g <- genes(txdb)
g[gene.name]

t <- transcriptsBy(txdb, "gene")
t[gene.name]
## btw v14prep and v15, the "tx_name" column was fixed

eg <- exonsBy(txdb, "gene")
eg[gene.name]
## btw v15 and v17, the "exon_name" column was fixed

et <- exonsBy(txdb, "tx", use.names=TRUE)
names(et)[grep(gene.name, names(et))]
et[paste0(gene.name, ".t01")]

f <- fiveUTRsByTranscript(txdb, use.names=TRUE)
names(f)[grep(gene.name, names(f))]
f[paste0(gene.name, ".t01")]
