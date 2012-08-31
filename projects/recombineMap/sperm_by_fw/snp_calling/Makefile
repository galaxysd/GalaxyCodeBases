all: pkusnp  spermsnp

spermsnp: call_sperm_snp.cpp gzstream.cpp
	g++ -O3 -o $@  $^ -lz

pkusnp: call_diploid_snp.cpp gzstream.cpp
	g++ -O3 -o $@  $^ -lz

clean:
	rm  pkusnp  spermsnp
