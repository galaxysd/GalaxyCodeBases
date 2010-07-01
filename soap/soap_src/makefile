CC=	g++

FLAGS= -O3 -DREAD_36 -static -DMAXGAP=3 -DMAXHITS=1000 -DDB_CHR -DTHREAD

ifdef DEBUG
FLAGS= -pg -O3 -DREAD_48 -static -DMAXGAP=3 -DMAXHITS=1000 -DDB_CHR -DTHREAD
endif

#-ggdb  -mpowerpc64 -maix64 #-mabi=64 # -qipa -q64
#DB_CHR  DB_CONTIG  DB_SHORT  DB_HUGE
#-DREAD_36 -DREAD_48

LIBS=  
THREAD=	-lpthread

SRC=	$(sort $(wildcard *.cpp))
OBJS=	$(subst .cpp, .o, $(SRC))

all: soap pes_dealign kmerfreq kmerfreq_vsref

%.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) -c $< -o $@

soap: main.o param.o align.o pairs.o dbseq.o reads.o utilities.o
	$(CC) $(FLAGS) $(LIBS) $^ -o $@ $(THREAD)

pes_dealign: pes_dealign.o dealign.o param.o
	$(CC) $(FLAGS) $(LIBS) $^ -o $@

kmerfreq: kmerfreq.o
	$(CC) $(FLAGS) $(LIBS) $^ -o $@
kmerfreq_vsref: kmerfreq_vsref.o
	$(CC) $(FLAGS) $(LIBS) $^ -o $@
	
clean:
	rm -f *.o soap pes_dealign kmerfreq kmerfreq_vsref
install:
	chmod 755 *.pl
	cp soap pes_dealign kmerfreq kmerfreq_vsref *.pl  ../ #/public/home/lirq/disk5/solexa/bin/
	
