DEFINE =
CXX = g++
CXXFLAGS = -fomit-frame-pointer -O3 -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0  #-MMD -MP -MF #-g3 -Wall -maccumulate-outgoing-args

LFLAGS = -lz

all: soapsnp
.PHONY: all

objects: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o
$(objects): %.o: soap_snp.h makefile

soapsnp: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o makefile
	$(CXX) $(CXXFLAGS) call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o -o soapsnp $(LFLAGS)

.PHONY: clean
clean:
	rm -f *.o soapsnp
