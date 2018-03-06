#bin/sh

g++ -O2 -ansi -pedantic -W -Wall -Wshadow -I../pythia8183/include pythielsenPART276.cc -o pythiatest.exe \
	-L../pythia8183/lib/archive -lpythia8 -llhapdfdummy
