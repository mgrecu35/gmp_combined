CC = gcc -Wall -g
CFLAGS = -O2
OBJECTS = main.o
INCLUDES = -I../include
LIBS = -L../lib -lmultiscatter \
	-lm -lgsl -lgslcblas
PROGRAM = test_blind_test

ms:
	$(CC) $(CFLAGS) $(INCLUDES) -c $(PROGRAM).c
	$(CC) $(CFLAGS) $(PROGRAM).o -o $(PROGRAM) $(LIBS)

test: test.c
	gcc $(INCLUDES) -c test.c
	gcc test.o -o test $(LIBS)