MEXEXT:=`mexext`
FC:=/usr/local/bin/gfortran

all : calc1ints.$(MEXEXT) calc2ints.$(MEXEXT)

calc1ints.$(MEXEXT) : calc1ints.f90 types.f90
	mex LD=g++ FC=$(FC) FFLAGS=-cpp -I./gen1int/include -I./gen1int/build/modules -L./gen1int/build -lgen1int -output calc1ints types.f90 calc1ints.f90

calc2ints.$(MEXEXT) : calc2ints.cc
	mex CC=c++ CXX=c++ LD=c++ -I./libint/include -L./libint/lib/.libs -lint2 calc2ints.cc

rebuild : 
	touch calc1ints.f90 calc2ints.cc
	make

clean :
	-rm calc1ints.$(MEXEXT)
	-rm calc2ints.$(MEXEXT)
	-rm types.mod

