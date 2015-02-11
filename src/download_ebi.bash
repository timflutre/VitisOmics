#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/vitis_vinifera

wget --timestamping ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/vitis_vinifera/dna/Vitis_vinifera.IGGP_12x.21.dna.genome.fa.gz

wget --timestamping ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/vitis_vinifera/dna/Vitis_vinifera.IGGP_12x.21.dna_sm.genome.fa.gz
