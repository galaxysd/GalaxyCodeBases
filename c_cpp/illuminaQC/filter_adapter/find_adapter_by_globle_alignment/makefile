PRGS = find_adapter_simple find_adapter
CC = g++
CFLAGS = -O3

all: $(PRGS)

find_adapter_simple: common.cpp find_adapter_simple.cpp locate_adapter.cpp gzstream.cpp
	$(CC) $(CFLAGS) -o $@   $^ -lz

find_adapter: common.cpp find_adapter.cpp locate_adapter.cpp gzstream.cpp
	$(CC) $(CFLAGS) -o  $@  $^ -lz

clean:
	rm $(PRGS)
