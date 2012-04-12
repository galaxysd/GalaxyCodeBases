DEFINE =

CXX = g++ #/ifs1/cloud/APT/wangbq/Devel/intel/Compiler/11.1/069/bin/intel64/icpc

DEBUG = -g

CXXFLAGS = -fomit-frame-pointer -O3 -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0 #-pipe -freg-struct-return   -maccumulate-outgoing-args -m64 -ftree-loop-linear -fprefetch-loop-arrays -fno-optional-diags -fno-threadsafe-statics -mfpmath=sse -m128bit-long-double #-MMD -MP -MF #-g3 -Wall -maccumulate-outgoing-args -Wstrict-aliasing=2 -fstack-check

DFLAGS= #-static -D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -D_CURSES_LIB=1

LIBPATH =  -I/ifs2/CD/tangzhb/bin/boost/include/ -L/panfs/CD/tangzhb/lib/SOAPsnpLib # -L/ifs2/CD/tangzhb/bin/zlib-1.2.5/

LFLAGS = -lz -lbam  -lm -lpthread #-liomp5 #-lcurses
# libiomp5 of Intel not used at all. Nor did -lgomp

objects = call_genotype.o chromosome.o matrix.o normal_dis.o prior.o rank_sum.o main.o tools.o SamCtrl.o MatrixManager.o Readwin.o FileListManager.o Call_winManager.o accessControl.o SfsMethod.o Files.o ThreadManager.o CThreadPool.o gzstream.o

all: DEBUG
.PHONY: all


$(objects): %.o:

DEBUG: CXXFLAGS +=  $(DFLAGS)
DEBUG: $(objects)
	$(CXX) $(CXXFLAGS) $(LIBPATH) $(objects) -o soapsnp $(LFLAGS)

release: CXX = /ifs1/cloud/APT/wangbq/Devel/intel/Compiler/11.1/069/bin/intel64/icpc
release: CXXFLAGS = -O3 -parallel  -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 $(DFLAGS) #-opt-report 1 -opt-report-phase=ipo 
release: $(objects)
	$(CXX) $(CXXFLAGS) $(LIBPATH) $(objects) -o soapsnp $(LFLAGS)

.PHONY: clean
clean:
	rm -f *.o soapsnp
