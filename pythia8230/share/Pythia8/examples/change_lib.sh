#!/bin/bash

folder=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/lib
pdf=$1

cd $folder

if [ $1 = "MED" ]; then
  echo "for MED"
  rm libpythia8lhapdf5.so
  cp med_libpythia8lhapdf5.so libpythia8lhapdf5.so
else
  echo "for VAC"
  rm libpythia8lhapdf5.so
  cp vac_libpythia8lhapdf5.so libpythia8lhapdf5.so
fi


