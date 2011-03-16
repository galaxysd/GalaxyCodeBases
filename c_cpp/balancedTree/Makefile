ARG = -Wall -O3 -pipe -march=core2 -mtune=generic -I. ./libpspp/bt.c

all: bin test
bin:
	gcc $(ARG) bt-test.c -o bt-test
	strip bt-test
	gcc $(ARG) bt-tmp.c -o btmp
	@echo "[!]Make done."
test:
	@echo "[!]Begin test ..."
	./bt-test insert-any-remove-any
	./bt-test insert-any-remove-same
	./bt-test insert-any-remove-reverse
	./bt-test random-sequence
	./bt-test insert-ordered
	./bt-test find-ge-le
	./bt-test moved
	./bt-test changed

d:
	gcc $(ARG) bt-tmp.c -o btmp

