fCSV = $(wildcard _TabIn/*.csv)
fTXT = $(wildcard _TabIn/*.txt)

report:
	./doAll.r -i $(fCSV) -m $(fTXT)
	@echo =================
	@cat pdf/_QC.txt
	@echo =================

clean:
	-rm *.log
	-rm -fr out

reset: clean
	-rm _TabIn/*
	-rm -fr pdf

pkg:
	tar -C .. -czvf ../rmd4.tgz rmd4/20200608rep.Rmd rmd4/doAll.r rmd4/files rmd4/inc rmd4/_TabIn rmd4/Makefile rmd4/_TabIn.example
