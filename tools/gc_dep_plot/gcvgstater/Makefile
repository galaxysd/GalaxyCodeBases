OBJ = gzstream.o self_util.o main.o stat_soap_coverage.o
CC = g++
LIBRARY = -O3 -lz -Wall -pipe -march=core2 -mtune=generic -pedantic -Wfloat-equal

gc_coverage_bias:${OBJ}
	${CC} -o gc_coverage_bias ${LIBRARY} ${OBJ}
gzstream.o:gzstream.h gzstream.cpp
	${CC} -c ${LIBRARY} gzstream.h gzstream.cpp
self_util.o:self_util.cpp self_util.h
	${CC} -c ${LIBRARY} self_util.cpp self_util.h
main.o:main.cpp
	${CC} -c ${LIBRARY} main.cpp
stat_soap_coverage.o:stat_soap_coverage.h stat_soap_coverage.cpp
	${CC} -c ${LIBRARY} stat_soap_coverage.h stat_soap_coverage.cpp
clean:
	-rm ${OBJ} *.gch
