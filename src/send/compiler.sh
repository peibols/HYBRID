#!/bin/sh

g++ -std=c++11 main.cc JetSafe.cc HadWake.cc Wake.cc dEdx.cc SonsMomenta.cc Eloss.cc Quench.cc Glauber.cc Hydro.cc Random.cc Parton.cc \
	 ../pythia8223/lib/libpythia8.a -o main -I../pythia8223/include \
	`../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` \
	-O2 -ansi -pedantic -W -Wall -Wshadow -Wcast-align -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winvalid-pch -Werror=missing-braces -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wswitch-default -Wno-unused -Wvariadic-macros -Wwrite-strings -Werror=declaration-after-statement -Werror=implicit-function-declaration -Werror=nested-externs -Werror=old-style-definition -Werror=strict-prototypes -fPIC -Wl,-rpath,../pythia8223/lib -ldl
