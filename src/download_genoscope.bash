#!/usr/bin/env bash

# Copyright (C) 2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/

mkdir -p 8X
cd 8X/

mkdir -p annotation
cd annotation/

if [ ! -f Vitis_vinifera_annotation_v1.gff.gz ]; then
  wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/8X/annotation/Vitis_vinifera_annotation_v1.gff
  gzip Vitis_vinifera_annotation_v1.gff
fi

cd .. # up from annotation/

cd .. # up from 8X/


mkdir -p 12X
cd 12X/

mkdir -p annotation
cd annotation/

if [ ! -f Vitis_vinifera_annotation.gff.gz ]; then
  wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/annotation/Vitis_vinifera_annotation.gff
  gzip Vitis_vinifera_annotation.gff
fi

cd .. # up from annotation/

cd .. # up from 12X/
