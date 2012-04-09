g++ -O3 -o filter_PEreads filter_adapter_lowqual_PE.cpp gzstream.cpp -lz
g++ -O3 -o filter_SEreads filter_adapter_lowqual_SE.cpp gzstream.cpp -lz
time  ./filter_raw_reads -q1 -a1 Ecoli_haplid_25Xreads_100_500_1.fq.gz Ecoli_haplid_25Xreads_100_500_2.fq.gz 1.clean 2.clean stat 

g++ -O3 -o calculate_duplcation duplication.cpp gzstream.cpp -lz

