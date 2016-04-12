#!/usr/bin/env bash

# Copyright (C) 2015-2016 Institut National de la Recherche Agronomique (INRA)
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

mkdir -p assembly
cd assembly/

wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/8X/assembly/Readme_First

cd .. # up from assembly/

cd .. # up from 8X/


mkdir -p 12X
cd 12X/

mkdir -p annotation
cd annotation/

wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/annotation/Vitis_vinifera_annotation.gff.gz

cd .. # up from annotation/

mkdir -p assembly
cd assembly/

mkdir -p goldenpath
cd goldenpath/
wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/Readme_First
wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/chr.agp
wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/chr.agp.info
wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/chr.lg
wget --timestamping http://www.genoscope.cns.fr/externe/Download/Projets/Projet_ML/data/12X/assembly/goldenpath/Markers.gff.gz
cd .. # from goldenpath/

cd .. # up from assembly/

cd .. # up from 12X/
