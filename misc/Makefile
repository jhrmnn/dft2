all: tprod

tprod:
	$(eval tmpfile := $(shell mktemp))
	wget -O ${tmpfile} http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/16275/versions/5/download/zip
	mkdir tprod
	unzip -d tprod ${tmpfile}
	rm ${tmpfile}

clean:

distclean: clean
	-rm -r tprod
