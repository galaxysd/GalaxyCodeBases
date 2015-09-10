#PROG=../bin/mvicuna
PROG = ../bin/alnmethly

#COMPILER=/opt/local/bin/g++
# server 
COMPILER=g++ 
FLAGs=-std=c++0x
OMP=-fopenmp

all:
	$(COMPILER) $(OMP) $(FLAGs) -O3 *.cpp -o $(PROG)
	
clean:
	rm -rf $(PROG)
