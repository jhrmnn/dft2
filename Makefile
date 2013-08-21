default :
	make deploy
	make install
	make test

unpack :
	tar xf deploy.tgz

remove :
	rm -r data
	rm dft/getLebedevSphere.m
	rm -r ints/gen1int
	rm -r ints/libint
	rm -r misc/tprod

make install :
	cd dft
	make
	cd ..
	cd ints
	cd genint
	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make
	cd ../..
	cd libint
	./configure
	make
	cd ..
	make
	cd ..

test :
	matlab -nojvm -r "test; exit"
