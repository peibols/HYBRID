#!/bin/bash

folder=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/include/Pythia8Plugins
pdf=$1

cd $folder

if [ $1 = "MED" ]; then
  echo "for MED"
  rm LHAPDF5.h
  cp my_LHAPDF5.h LHAPDF5.h
  cd ..
  cd ..
  make  
else
  echo "for VAC"
  rm LHAPDF5.h
  cp orig_LHAPDF5.h LHAPDF5.h
  cd ..
  cd ..
  make
fi


