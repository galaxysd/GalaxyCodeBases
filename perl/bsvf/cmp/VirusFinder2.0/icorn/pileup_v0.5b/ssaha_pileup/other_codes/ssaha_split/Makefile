CC = cc
CFLAGS = -O4

LIB = -lm

OBJS=	\
		ssaha_split.o \
                fast.o \

SRCS= $(OBJS:.o=.c)

INCS=   \
                fasta.h \

PROGRAM = ssaha_split

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
