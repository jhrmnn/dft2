mexext := $(shell mexext)

all: calc1ints.${mexext} calc2ints.${mexext}

calc1ints.${mexext}: types.f90 calc1ints.f90 | gen1int/build/libgen1int.a
	mex LD=g++ FC=gfortran FFLAGS=-cpp -I./gen1int/include -I./gen1int/build/modules -L./gen1int/build -lgen1int -output $@ $^

gen1int/build/libgen1int.a:
	mkdir -p gen1int && wget -O - https://gitlab.com/bingao/gen1int/repository/archive.tar.gz?ref=70d10a37390f4e42a90b2951cb3fec3e9a870d13 | tar -C gen1int -zx --strip-components 1
	mkdir -p $(dir $@)
	cd $(dir $@) && cmake -DCMAKE_BUILD_TYPE=Release ..
	${MAKE} -C $(dir $@)

calc2ints.${mexext}: calc2ints.cc | libint/lib/.libs/libint2.a
	mex LD=g++ CC=g++ -I./libint/include -I./libint/include/libint2 -L./libint/lib/.libs -lint2 $^

libint/lib/.libs/libint2.a:
	mkdir -p libint && wget -O - https://github.com/evaleev/libint/releases/download/v2.1.0-beta2/libint-2.1.0-beta2.tgz | tar -C libint -zx --strip-components 1
	cd libint && ./configure
	${MAKE} -C libint

clean:
	-rm *.${mexext}
	-rm *.mod

distclean: clean
	-rm -r gen1int/build
	-rm aclocal.m4
	-rm -r autom4te.cache
	-rm -r build
	-rm configure
