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

if ! hash makeblastdb 2>/dev/null; then
  echo "ERROR: can't find makeblastdb" 1>&2
  exit 1
fi

inFaGzFile=$1

if [ ! -f "${inFaGzFile}" ]; then
  echo "ERROR: can't find file ${inFaGzFile}" 1>&2
  exit 1
fi

baseName=$(basename ${inFaGzFile} .fa.gz)

zcat ${inFaGzFile} \
  | makeblastdb \
      -dbtype nucl \
      -in - \
      -title ${baseName} \
      -out ${baseName} \
      -logfile makeblastdb_${baseName}.log

echo "done!"
date
