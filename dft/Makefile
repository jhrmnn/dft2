mexext := $(shell mexext)

all: eval_xc.${mexext} getLebedevSphere.m

eval_xc.${mexext}: xc_func.f eval_xc.f90
	mex LD=g++ FC=gfortran FFLAGS=-cpp -I../ints -output eval_xc xc_func.f eval_xc.f90 

getLebedevSphere.m:
	$(eval tmpfile := $(shell mktemp))
	$(eval tmpdir := $(shell mktemp -d))
	wget -O ${tmpfile} http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/27097/versions/1/download/zip/getLebedevSphere.zip
	unzip -d ${tmpdir} ${tmpfile}
	mv ${tmpdir}/$@ .
	rm ${tmpfile}
	rm -r ${tmpdir}

xc_func.f: dfauto funcs/*.df
	cp funcs/*.df .
	FC=gfortran TMPDIR= sh dfauto -s -t a4paper -x 12pt -o xc_func.f *.df
	rm *.df
	dvipdfm -o funcs.pdf acg.dvi >/dev/null
	rm acg.*
	rm *.tex

dfauto: dfauto.patch
	wget -O - "http://cpc.cs.qub.ac.uk/cpc/cgi-bin/processrequest.pl/?catid=adny_v1_0&usertype=fromsummary&deliverytype=www&email=" | gunzip | tail -n +4 >$@
	patch $@ $<

clean: 
	cp funcs/*.df . 
	-FC=gfortran TMPDIR= sh dfauto -r -o xc_func.f *.df
	-rm *.df
	-rm eval_xc.${mexext}
	-rm funcs.pdf

distclean: clean
	-rm getLebedevSphere.m
	-rm dfauto
