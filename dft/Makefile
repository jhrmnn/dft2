MEXEXT:=`mexext`
FC:=/usr/local/bin/gfortran

all : eval_xc.$(MEXEXT)

xc_func.f : dfauto funcs/*.df
	cp funcs/*.df .
	FC=$(FC) TMPDIR= sh dfauto -s -t a4paper -x 12pt -o xc_func.f *.df
	rm *.df
	dvipdfm -o funcs.pdf acg.dvi > /dev/null
	rm acg.*
	mv *.tex funcs

eval_xc.$(MEXEXT) : xc_func.f eval_xc.f90
	mex LD=g++ FC=$(FC) FFLAGS=-cpp -I../ints -output eval_xc xc_func.f eval_xc.f90 

clean : 
	cp funcs/*.df .
	FC=$(FC) TMPDIR= sh dfauto -r -o xc_func.f *.df
	rm *.df

rebuild :
	touch dfauto
	make
