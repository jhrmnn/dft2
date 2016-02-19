all:
	make -C ints
	make -C dft
	make -C misc

clean:
	make -C ints clean
	make -C dft clean
	make -C misc clean

test:
	matlab -nojvm -r "test; exit"
