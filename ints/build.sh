mex LD=g++ FC=/usr/local/bin/gfortran FFLAGS=-cpp -I./gen1int/include -I./gen1int/build/modules -L./gen1int/build -lgen1int -output calc1ints types.f90 calc1ints.f90

mex CC=c++ CXX=c++ LD=c++ -I./libint-2.0.3-stable/include -L./libint-2.0.3-stable/lib/.libs -lint2 calc2ints.cc

