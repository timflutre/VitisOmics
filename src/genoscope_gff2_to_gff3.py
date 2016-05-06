#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: convert GFF2 files from Genoscope for V. vinifera to GFF3
# Copyright (C) 2016 Timothée Flutre
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: https://github.com/timflutre/VitisOmics

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import sys
import os
import getopt
import time
import datetime
from subprocess import Popen, PIPE
import math
import gzip
import re

from Bio import SeqIO

if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        msg = "ERROR: Python should be in version 2.7 or higher"
        sys.stderr.write("%s\n\n" % msg)
        sys.exit(1)
        
progVersion = "1.0.0" # http://semver.org/


# http://stackoverflow.com/a/3033342/597069
def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


class GenoscopeGff2ToGff3(object):
    
    def __init__(self):
        self.verbose = 1
        self.gff2File = None
        self.faFile = "VITVI_PN40024_12x_v0_chroms_URGI.fa.gz"
        self.gff3File = None
        
        
    def help(self):
        """
        Display the help on stdout.
        
        The format complies with help2man (http://www.gnu.org/s/help2man)
        """
        msg = "`%s' converts GFF2 files from Genoscope for V. vinifera to GFF3.\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Usage: %s [OPTIONS] ...\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Options:\n"
        msg += "  -h, --help\tdisplay the help and exit\n"
        msg += "  -V, --version\toutput version information and exit\n"
        msg += "  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
        msg += "      --gff2\tpath to the input file in GFF2 (can be gzipped)\n"
        msg += "      --fa\tpath to the reference genome in fasta format (can be gzipped)\n"
        msg += "\t\tdefault=VITVI_PN40024_12x_v0_chroms_URGI.fa.gz\n"
        msg += "      --gff3\tpath to the output file in GFF3 (can be gzipped)\n"
        msg += "\n"
        msg += "Examples:\n"
        msg += "  %s --gff2 Vitis_vinifera_annotation.gff.gz --gff3 Vitis_vinifera_annotation.gff3.gz\n" % os.path.basename(sys.argv[0])
        msg += "\n"
        msg += "Report bugs to <timothee.flutre@supagro.inra.fr>."
        print(msg); sys.stdout.flush()
        
        
    def version(self):
        """
        Display version and license information on stdout.
        
        The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012).
        """
        msg = "%s %s\n" % (os.path.basename(sys.argv[0]), progVersion)
        msg += "\n"
        msg += "Copyright (C) 2016 Timothée Flutre.\n"
        msg += "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
        msg += "\n"
        msg += "Written by Timothée Flutre [cre,aut]."
        print(msg.encode("utf8")); sys.stdout.flush()
        
        
    def setAttributesFromCmdLine(self):
        """
        Parse the command-line arguments.
        """
        try:
            opts, args = getopt.getopt( sys.argv[1:], "hVv:",
                                        ["help", "version", "verbose=",
                                         "gff2=", "fa=", "gff3="])
        except getopt.GetoptError as err:
            sys.stderr.write("%s\n\n" % str(err))
            self.help()
            sys.exit(2)
        for o, a in opts:
            if o == "-h" or o == "--help":
                self.help()
                sys.exit(0)
            elif o == "-V" or o == "--version":
                self.version()
                sys.exit(0)
            elif o == "-v" or o == "--verbose":
                self.verbose = int(a)
            elif o == "--gff2":
                 self.gff2File = a
            elif o == "--fa":
                self.faFile = a
            elif o == "--gff3":
                 self.gff3File = a
            else:
                assert False, "invalid option"
                
                
    def checkAttributes(self):
        """
        Check the values of the command-line parameters.
        """
        if not self.gff2File:
            msg = "ERROR: missing compulsory option --gff2"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.gff2File):
            msg = "ERROR: can't find file %s" % self.gff2File
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.faFile:
            msg = "ERROR: missing compulsory option --fa"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not os.path.exists(self.faFile):
            msg = "ERROR: can't find file %s" % self.faFile
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
        if not self.gff3File:
            msg = "ERROR: missing compulsory option --gff3"
            sys.stderr.write("%s\n\n" % msg)
            self.help()
            sys.exit(1)
            
            
    def writeHeader(self, gff2Handle, faHandle, gff3Handle):
        if self.verbose > 0:
            print("write header ...")
            sys.stdout.flush()
            
        txt = "##gff-version 3"
        gff3Handle.write("%s\n" % txt)
        
        ## load chromosomes from fasta file
        dChrsFa = SeqIO.to_dict(SeqIO.parse(faHandle,"fasta"))
        
        ## load chromosomes from gff2 file
        dChrsGff2 = {}
        while True:
            line = gff2Handle.readline()
            if line == "":
                break
            tokens = line.split("\t")
            if tokens[0] not in dChrsGff2:
                dChrsGff2[tokens[0]] = [1, 1]
            dChrsGff2[tokens[0]][1] = max(dChrsGff2[tokens[0]][1], int(tokens[4]))
            
        ## check all chromosomes are present in both files
        for chrom in dChrsFa:
            if chrom not in dChrsGff2:
                msg = "chromosome %s from fasta file absent from gff2 file" % chrom
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
        for chrom in dChrsGff2:
            if chrom not in dChrsFa:
                msg = "chromosome %s from gff2 file absent from fa file" % chrom
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
                
        ## check all annotations have coordinates inside chrom lengths
        for chrom in dChrsGff2:
            if dChrsGff2[chrom][1] > len(dChrsFa[chrom]):
                msg = "an annotation on %s has coordinate larger than %i" \
                      % len(dChrsFa[chrom])
                sys.stderr.write("%s\n" % msg)
                sys.exit(1)
                
        ## make gff3 header
        lChrs = dChrsFa.keys()
        lChrs.sort(key=natural_key)
        for chrom in lChrs:
            txt = "##sequence-region %s 1 %i" % (chrom, len(dChrsFa[chrom]))
            gff3Handle.write("%s\n" % txt)
            
        if self.verbose > 0:
            print("%i seqids" % len(dChrsFa))
            sys.stdout.flush()
            
            
    def convertAnnot(self, gff2Handle, gff3Handle):
        """
        Skip UTRs.
        """
        if self.verbose > 0:
            print("convert annotations ...")
            sys.stdout.flush()
            
        gff2Handle.seek(0)
        fmt = "\t".join(["{0}","{1}","{2}","{3}","{4}","{5}","{6}","{7}"])
        dGenes = {} # gene_ID: [number of mRNA(s), number of CDS(s)]
        
        while True:
            line = gff2Handle.readline()
            if line == "":
                break
            tokens = line.split("\t")
            tokensAttr = tokens[8].split(" ")
            
            txt = fmt.format(*(tokens[0:8])) # no change to the first 8 fields
            
            if tokens[2] == "gene":
                txt = "%s\tID=%s;Name=%s" % (txt, tokensAttr[1], tokensAttr[1])
                dGenes[tokensAttr[1]] = {"mRNA":0, "CDS":0}
                
            elif tokens[2] == "mRNA":
                geneId = tokensAttr[1].replace("GSVIVT", "GSVIVG")
                dGenes[geneId]["mRNA"] += 1
                txId = "%s" % (tokensAttr[1])
                txt = "%s" % txt
                txt += "\tID=%s" % txId
                txt += ";Name=%s" % txId
                txt += ";Parent=%s" % geneId
                
            elif tokens[2] == "CDS":
                geneId = tokensAttr[1].replace("GSVIVT", "GSVIVG")
                dGenes[geneId]["CDS"] += 1
                txId = tokensAttr[1]
                cdsId = "%s.cds.%d" % (txId, dGenes[geneId]["CDS"])
                txt = "%s" % txt
                txt += "\tID=%s" % cdsId
                txt += ";Name=%s" % cdsId
                txt += ";Parent=%s" % txId
                
            gff3Handle.write("%s\n" % txt)
            
        if self.verbose > 0:
            print("%i genes" % len(dGenes))
            nbTx = 0
            for gene in dGenes:
                nbTx += dGenes[gene]["mRNA"]
            print("%i transcripts" % nbTx)
            nbCds = 0
            for gene in dGenes:
                nbCds += dGenes[gene]["CDS"]
            print("%i transcripts" % nbCds)
            sys.stdout.flush()
            
            
    def run(self):
        if self.gff2File.endswith(".gz"):
            gff2Handle = gzip.open(self.gff2File, "r")
        else:
            gff2Handle = open(self.gff2File, "r")
        if self.faFile.endswith(".gz"):
            faHandle = gzip.open(self.faFile, "r")
        else:
            faHandle = open(self.faFile, "r")
        if self.gff3File.endswith(".gz"):
            gff3Handle = gzip.open(self.gff3File, "w")
        else:
            gff3Handle = open(self.gff3File, "w")
            
        self.writeHeader(gff2Handle, faHandle, gff3Handle)
        self.convertAnnot(gff2Handle, gff3Handle)
        
        gff2Handle.close()
        faHandle.close()
        gff3Handle.close()
        
        
if __name__ == "__main__":
    i = GenoscopeGff2ToGff3()
    
    i.setAttributesFromCmdLine()
    
    i.checkAttributes()
    
    if i.verbose > 0:
        startTime = time.time()
        msg = "START %s %s %s" % (os.path.basename(sys.argv[0]),
                                  progVersion,
                                  time.strftime("%Y-%m-%d %H:%M:%S"))
        msg += "\ncmd-line: %s" % ' '.join(sys.argv)
        msg += "\ncwd: %s" % os.getcwd()
        print(msg); sys.stdout.flush()
        
    i.run()
    
    if i.verbose > 0:
        msg = "END %s %s %s" % (os.path.basename(sys.argv[0]),
                                progVersion,
                                time.strftime("%Y-%m-%d %H:%M:%S"))
        endTime = time.time()
        runLength = datetime.timedelta(seconds=
                                       math.floor(endTime - startTime))
        msg += " (%s" % str(runLength)
        if "linux" in sys.platform:
            p = Popen(["grep", "VmHWM", "/proc/%s/status" % os.getpid()],
                      shell=False, stdout=PIPE).communicate()
            maxMem = p[0].split()[1]
            msg += "; %s kB)" % maxMem
        else:
            msg += ")"
        print(msg); sys.stdout.flush()
