#!/bin/bash

PYTHIA_BIN=/home/peibols/projects/rrg-jeon-ac/peibols/software/pythia-install/bin
PYTHIA_INCLUDE=/home/peibols/projects/rrg-jeon-ac/peibols/software/pythia-install/include
PYTHIA_LIB=/home/peibols/projects/rrg-jeon-ac/peibols/software/pythia-install/lib
PREFIX_SHARE=/home/peibols/projects/rrg-jeon-ac/peibols/software/pythia-install/share/Pythia8

FASTJET3_BIN=/home/peibols/projects/rrg-jeon-ac/peibols/software/fastjet-install/bin
FASTJET3_INCLUDE=/home/peibols/projects/rrg-jeon-ac/peibols/software/fastjet-install/include
FASTJET3_LIB=/home/peibols/projects/rrg-jeon-ac/peibols/software/fastjet-install/lib

gcc -g -std=c++0x -mcmodel=medium \
	$1.cc  shower_analysis.cc Tree.cc JetObs.cc JetSafe.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber_$2.cc Hydro_$2.cc Random.cc Parton.cc Hadron.cc Lund.cc -o $1_$2 \
	${PYTHIA_LIB}/libpythia8.a -I${PYTHIA_INCLUDE} -L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl \
	-I${FASTJET3_INCLUDE} -L${FASTJET3_LIB} -Wl,-rpath,${FASTJET3_LIB} -lfastjet `${FASTJET3_BIN}/fastjet-config --cxxflags --libs --plugins` \
	-O2 -ansi -pedantic -W -Wall -Wshadow -Wcast-align -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winvalid-pch -Werror=missing-braces -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wno-unused -Wvariadic-macros -Wwrite-strings -Werror=declaration-after-statement -Werror=implicit-function-declaration -Werror=nested-externs -Werror=old-style-definition -Werror=strict-prototypes -fPIC
