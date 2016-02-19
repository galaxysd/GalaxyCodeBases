#PROG=../bin/mvicuna
PROG = ../bin/alnmethly

#COMPILER=/opt/local/bin/g++
# server 
COMPILER=g++-5 
FLAGs=-std=c++0x
OMP=-fopenmp

all:
	mkdir ../bin
	$(COMPILER) $(OMP) $(FLAGs) -O3 *.cpp -o $(PROG)
	ln -s ../t.pl ../bin
	
clean:
	rm -rf $(PROG)
	rm -fr ../bin/
