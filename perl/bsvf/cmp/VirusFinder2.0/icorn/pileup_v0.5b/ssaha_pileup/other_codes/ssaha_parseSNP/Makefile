CC = gcc
CFLAGS = -Wall -O4

LIB = -lm

OBJS=	\
                fast.o \
		ssaha_parseSNP.o \

SRCS= $(OBJS:.o=.c)

PROGRAM = ssaha_parseSNP 

$(PROGRAM): $(OBJS) $(INCS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJS) $(LIB)

clean:
	rm -f $(PROGRAM) *.o
