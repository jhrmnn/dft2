default :
	make unpack
	make install
	make test

unpack :
	tar xf deploy.tgz

clean :
	-rm -r data
	-rm dft/getLebedevSphere.m
	make -C dft clean
	-rm -r ints/gen1int
	-rm -r ints/libint
	make -C ints clean
	-rm -r misc/tprod

make install :
	mkdir ints/gen1int/build
	cd ints/gen1int/build; cmake -DCMAKE_BUILD_TYPE=Release ..
	make -C ints/gen1int/build
	cd ints/libint; ./configure
	make -C ints/libint
	make -C ints
	make -C dft rebuild

test :
	matlab -nojvm -r "test; exit"
