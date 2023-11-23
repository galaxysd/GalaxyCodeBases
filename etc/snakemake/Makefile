run: clean
	snakemake -p -c4 >r.log 2>r.err
	@echo ... ...
	cat -n r.log
	@echo ... ...
	cat -n r.err | head -n26
	@echo ... ...
	cat -n r.err | tail -n19
	@echo ... ...
	ls -l results.sha2 r.log r.err
	@echo ... ...
	find .Trash/

dry:
	snakemake -p -c4 -n

help:
	snakemake -c2 help

clean:
	rm -fr r.log r.err .Trash results

pkg: clean
	rm -fr snakemakeExample.tgz
	tar -czvf snakemakeExample.tgz config workflow exampleData files.lst Makefile
	ls -l snakemakeExample.tgz
