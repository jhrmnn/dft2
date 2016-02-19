mexext := $(shell mexext)

all: calc1ints.${mexext} calc2ints.${mexext}

calc1ints.${mexext}: types.f90 calc1ints.f90 | gen1int/build/libgen1int.a
	mex LD=g++ FC=gfortran FFLAGS=-cpp -I./gen1int/include -I./gen1int/build/modules -L./gen1int/build -lgen1int -output $@ $^

gen1int/build/libgen1int.a:
	mkdir -p $(dir $@)
	cd $(dir $@) && cmake -DCMAKE_BUILD_TYPE=Release ..
	${MAKE} -C $(dir $@)

calc2ints.${mexext}: calc2ints.cc | libint/build/lib/.libs/libint2.a
	mex LD=g++ CC=g++ -I./libint/include -I./libint/include/libint2 -I./libint/build/include -I./libint/build/include/libint2 -L./libint/build/lib/.libs -lint2 $^

libint/lib/.libs/libint2.a:
	cd libint && ./autogen.sh && mkdir -p build && cd build && ../configure
	${MAKE} -C libint/build

clean:
	-rm *.${mexext}
	-rm *.mod

distclean: clean
	-rm -r gen1int/build
	-rm aclocal.m4
	-rm -r autom4te.cache
	-rm -r build
	-rm configure
