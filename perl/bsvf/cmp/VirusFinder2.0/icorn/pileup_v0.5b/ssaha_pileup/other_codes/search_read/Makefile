CC = gcc
CFLAGS = -Wall -O4

LIB = -lm

OBJS=	\
		search_read.o \
                fast.o \

SRCS= $(OBJS:.o=.c)

INCS=   \
                fasta.h \

PROGRAM = search_read

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
