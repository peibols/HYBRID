#!/bin/bash

export LD_LIBRARY_PATH=/software/CentOS-6/libraries/boost-1.57/lib:/software/CentOS-6/libraries/GSL/1.15/lib:/software/CentOS-6/tools/openmpi-1.8.3-intel/lib:/opt/torque/x86_64/lib:/software/CentOS-6/libraries/OpenBLAS_LAPACK/0.2.12-openmp-gcc/lib:/software/CentOS-6/tools/python-2.7.3/lib:/software/CentOS-6/tools/wx-2.8.12/lib:/software/CentOS-6/tools/hdf5/1.8.13-intel/lib:/software/CentOS-6/compilers/gcc-4.9.1/lib64:/software/CentOS-6/compilers/gcc-4.9.1/lib:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/compiler/lib/intel64:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl/lib/intel64:/gs/project/cqn-654-ad/peibols/lhapdf_pythia/lib

N=$1
pthatmin=1
seed=0

./change_lib.sh VAC
./vac_pdf_test $pthatmin $N $seed VAC
echo "Done VAC"

./change_lib.sh MED
./med_pdf_test $pthatmin $N $seed MED
echo "Done MED"
