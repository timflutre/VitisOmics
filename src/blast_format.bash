#!/usr/bin/env bash

# Copyright (C) 2014-2015 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Author: Timothée Flutre

# http://www.ncbi.nlm.nih.gov/books/NBK1762/
# http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/doc/blast/formatdb.html

set -e -o pipefail -u # http://stackoverflow.com/a/69808/597069

date
echo "start ..."

if [ $# -eq 0 ]; then
  echo "ERROR: missing input file (.fa.gz)" 1>&2
  exit 1
fi

inFaGzFile=$1

if [ ! -f "${inFaGzFile}" ]; then
  echo "ERROR: can't find file ${inFaGzFile}" 1>&2
  exit 1
fi

baseName=$(basename ${inFaGzFile} .fa.gz)

zcat ${inFaGzFile} \
  | formatdb \
      -i stdin \
      -l formatdb_${baseName}.log \
      -p F \
      -n ${baseName}

echo "done!"
date
