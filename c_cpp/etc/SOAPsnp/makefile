DEFINE =
CXX = g++
CXXFLAGS = -g -fomit-frame-pointer -O3 -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0  -I/ifs2/CD/tangzhb/bin/zlib-1.2.5/samtool/ -L/ifs2/CD/tangzhb/bin/zlib-1.2.5/#-MMD -MP -MF #-g3 -Wall -maccumulate-outgoing-args
#LIBPATH = -I/panfs/RD/tangzhb/bin/zlib-1.2.5/samtool/ -L/panfs/RD/tangzhb/bin/zlib-1.2.5/

LFLAGS = -lz -lbam -lcurses -lm

all: soapsnp
.PHONY: all

objects: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o tools.o
$(objects): %.o: soap_snp.h makefile

soapsnp: call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o tools.o SamCtrl.o makefile
	$(CXX) $(CXXFLAGS) $(LIBPATH) call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o tools.o SamCtrl.o -o soapsnp $(LFLAGS)

.PHONY: clean
clean:
	rm -f *.o soapsnp
