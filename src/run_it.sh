#!/bin/bash

#Args: VAC or MED pdf, pPb or PbPb system

nhat=5
N=30000
cent="00_05"
kappa=0.
alpha=0.442
tmethod=1

pdf=$1
system=$2

local="$PWD"
folder=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/lib

cd $folder

#SELECT PDF
if [ $1 = "MED" ] && [ $2 = "AuAu" ]; then
  echo "for AuAu MED"
  #rm libpythia8lhapdf5.so
  #cp Au_libpythia8lhapdf5.so libpythia8lhapdf5.so
elif [ $1 = "MED" ] && [ $2 = "pPb" ]; then
  echo "for pPb MED"
  #rm libpythia8lhapdf5.so
  #cp Pb_libpythia8lhapdf5.so libpythia8lhapdf5.so
elif [ $1 = "MED" ] && [ $2 = "PbPb" ]; then
  echo "for PbPb MED"
  #rm libpythia8lhapdf5.so
  #cp Pb_libpythia8lhapdf5.so libpythia8lhapdf5.so
else
  echo "for VAC"
  #rm libpythia8lhapdf5.so
  #cp p_libpythia8lhapdf5.so libpythia8lhapdf5.so
fi

cd $local

export LD_LIBRARY_PATH=/software/CentOS-6/libraries/boost-1.57/lib:/software/CentOS-6/libraries/GSL/1.15/lib:/software/CentOS-6/tools/openmpi-1.8.3-intel/lib:/opt/torque/x86_64/lib:/software/CentOS-6/libraries/OpenBLAS_LAPACK/0.2.12-openmp-gcc/lib:/software/CentOS-6/tools/python-2.7.3/lib:/software/CentOS-6/tools/wx-2.8.12/lib:/software/CentOS-6/tools/hdf5/1.8.13-intel/lib:/software/CentOS-6/compilers/gcc-4.9.1/lib64:/software/CentOS-6/compilers/gcc-4.9.1/lib:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/compiler/lib/intel64:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl/lib/intel64:/gs/project/cqn-654-ad/peibols/lhapdf_pythia/lib:/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/lib

export PYTHIA8DATA=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/share/Pythia8/xmldoc

export LHAPATH=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/share/lhapdf/PDFsets

./main_$2 $nhat $N $cent $kappa $alpha $tmethod $system $pdf
