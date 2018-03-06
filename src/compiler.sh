#!/bin/sh

system=$2

if [ "$2" == "PbPb" ]; then
	g++ -g -std=c++0x -mcmodel=medium $1.cc JetObs.cc JetSafe.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber_AA.cc Hydro_AA.cc Random.cc Parton.cc Hadron.cc Lund.cc \
		../../lhapdf_pythia/pythia8230/lib/libpythia8.a -o $1_PbPb -I../../lhapdf_pythia/pythia8230/include \
		`../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` \
		-O2 -ansi -pedantic -W -Wall -Wshadow -Wcast-align -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winvalid-pch -Werror=missing-braces -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wno-unused -Wvariadic-macros -Wwrite-strings -Werror=declaration-after-statement -Werror=implicit-function-declaration -Werror=nested-externs -Werror=old-style-definition -Werror=strict-prototypes -fPIC -Wl,-rpath,../pythia8223/lib -ldl
fi

if [ "$2" == "pPb" ]; then
	g++ -g -std=c++0x -mcmodel=medium $1.cc JetObs.cc JetSafe.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber_pA.cc Hydro_pA.cc Random.cc Parton.cc Hadron.cc Lund.cc \
                ../../lhapdf_pythia/pythia8230/lib/libpythia8.a -o $1_pPb -I../../lhapdf_pythia/pythia8230/include \
		`../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` \
                -O2 -ansi -pedantic -W -Wall -Wshadow -Wcast-align -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winvalid-pch -Werror=missing-braces -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wno-unused -Wvariadic-macros -Wwrite-strings -Werror=declaration-after-statement -Werror=implicit-function-declaration -Werror=nested-externs -Werror=old-style-definition -Werror=strict-prototypes -fPIC -Wl,-rpath,../pythia8223/lib -ldl
fi
