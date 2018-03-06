# PYTHIA configuration file.
# Generated on Tue Mar 6 16:58:27 EST 2018 with the user supplied options:
# --with-lhapdf5=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/

# Install directory prefixes.
PREFIX_BIN=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/bin
PREFIX_INCLUDE=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/include
PREFIX_LIB=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/lib
PREFIX_SHARE=/gs/project/cqn-654-ad/peibols/lhapdf_pythia/pythia8230/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=/software/CentOS-6/compilers/gcc-4.9.1/bin/g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=./
FASTJET3_LIB=./

# HEPMC2 configuration.
HEPMC2_USE=false
HEPMC2_BIN=
HEPMC2_INCLUDE=./
HEPMC2_LIB=./

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=true
LHAPDF5_BIN=/gs/project/cqn-654-ad/peibols/lhapdf_pythia//bin/
LHAPDF5_INCLUDE=/gs/project/cqn-654-ad/peibols/lhapdf_pythia//include
LHAPDF5_LIB=/gs/project/cqn-654-ad/peibols/lhapdf_pythia//lib

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=./
LHAPDF6_LIB=./

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=./
ROOT_LIB=./

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
