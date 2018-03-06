#!/bin/bash

#Args: VAC or MED pdf, pPb or PbPb system

#Nhat not meaningful right now
nhat=5

kappa=0.
alpha=0.43
tmethod=1

pdf=$1
system=$2

export LD_LIBRARY_PATH=/software/CentOS-6/libraries/boost-1.57/lib:/software/CentOS-6/libraries/GSL/1.15/lib:/software/CentOS-6/tools/openmpi-1.8.3-intel/lib:/opt/torque/x86_64/lib:/software/CentOS-6/libraries/OpenBLAS_LAPACK/0.2.12-openmp-gcc/lib:/software/CentOS-6/tools/python-2.7.3/lib:/software/CentOS-6/tools/wx-2.8.12/lib:/software/CentOS-6/tools/hdf5/1.8.13-intel/lib:/software/CentOS-6/compilers/gcc-4.9.1/lib64:/software/CentOS-6/compilers/gcc-4.9.1/lib:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/compiler/lib/intel64:/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl/lib/intel64:/gs/project/cqn-654-ad/peibols/lhapdf_pythia/lib:/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/lib

export PYTHIA8DATA=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/share/Pythia8/xmldoc

export LHAPATH=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/share/lhapdf/PDFsets

if [ $system = "PbPb" ]; then
	N=150000
	#0-5%
	cent="00_05"
	./main_$2 $nhat $N $cent $kappa $alpha $tmethod $system $pdf
	#0-5%
        cent="05_10"
        ./main_$2 $nhat $N $cent $kappa $alpha $tmethod $system $pdf

elif [ $system = "pPb" ]; then
	N=300000
        #0-10%
        cent="0-10"
        ./main_$2 $nhat $N $cent $kappa $alpha $tmethod $system $pdf
fi
