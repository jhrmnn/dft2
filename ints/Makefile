all : calc1ints.mexmaci64 calc2ints.mexmaci64

calc1ints.mex% : calc1ints.f90 types.f90
	mex LD=g++ FC=/usr/local/bin/gfortran FFLAGS=-cpp -I./gen1int/include -I./gen1int/build/modules -L./gen1int/build -lgen1int -output calc1ints types.f90 calc1ints.f90

calc2ints.mex% : calc2ints.cc
	mex CC=c++ CXX=c++ LD=c++ -I./libint-2.0.3-stable/include -L./libint-2.0.3-stable/lib/.libs -lint2 calc2ints.cc

rebuild : 
	touch calc1ints.f90 calc2ints.cc
	make

