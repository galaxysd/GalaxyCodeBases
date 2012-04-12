DEFINE =

CXX = g++

DEBUG = -g

CXXFLAGS = -Wstrict-aliasing=2 -O3 -fomit-frame-pointer -O3 -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0  #-MMD -MP -MF #-g3 -Wall -maccumulate-outgoing-args

LIBPATH = -I/ifs2/CD/luoruibang/bin/boost/include/  -L/ifs2/CD/tangzhb/bin/zlib-1.2.5/

LFLAGS = -lz -lbam -lcurses -lm -lpthread

objects = call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o tools.o SamCtrl.o MatrixManager.o Readwin.o FileListManager.o Call_winManager.o accessControl.o SfsMethod.o Files.o ThreadManager.o CThreadPool.o

all: soapsnp
.PHONY: all


$(objects): %.o:

soapsnp: $(objects)
	$(CXX) $(CXXFLAGS) $(LIBPATH) $(objects) -o soapsnp $(LFLAGS)

release: $(objects)
	$(CXX) $(CXXFLAGS) $(LIBPATH) $(objects) -o soapsnp $(LFLAGS)

.PHONY: clean
clean:
	rm -f *.o soapsnp
