ALL: main
.PHONY: ALL clean libraptor tests

CC = g++
LDFLAGS ?=
CFLAGS += -fPIC -O2

LIB_OBJS = Decoder.o Encoder.o Generators.o Helper.o StdAfx.o Symbol.o

Decoder.o: Decoder.cpp Decoder.h
	$(CC)  Decoder.cpp -o Decoder.o -c $(CFLAGS)
Encoder.o: Encoder.cpp Encoder.h
	$(CC) Encoder.cpp -o Encoder.o -c $(CFLAGS)
Generators.o: Generators.cpp Generators.h
	$(CC) Generators.cpp -o Generators.o -c $(CFLAGS)
Helper.o: Helper.cpp Helper.h
	$(CC) Helper.cpp -o Helper.o -c $(CFLAGS)
Main.o: Main.cpp
	$(CC) Main.cpp -o Main.o -c $(CFLAGS)
StdAfx.o: StdAfx.cpp StdAfx.h
	$(CC) StdAfx.cpp -o StdAfx.o -c $(CFLAGS)
Symbol.o: Symbol.cpp Symbol.h
	$(CC) Symbol.cpp -o Symbol.o -c $(CFLAGS)
test.o: test.cpp
	$(CC) test.cpp -o test.o -c $(CFLAGS)

main: $(LIB_OBJS) Main.o
	$(CC) $(LIB_OBJS) Main.o -o main

test_runner: $(LIB_OBJS) test.o
	$(CC) $(LIB_OBJS) test.o -o test_runner

tests: test_runner
	./test_runner

libraptorq: libraptorq.a
libraptorq.a: Decoder.o Encoder.o Generators.o Helper.o Symbol.o
	ar rcs libraptorq.a $^
#	$(CC) -shared -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o main test_runner libraptorq.so libraptorq.a libraptor.dll
