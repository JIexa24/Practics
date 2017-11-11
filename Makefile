CC = gcc
CFLAGS = -masm=intel -pthread

.PHONY: all COMPILE

all: COMPILE

COMPILE:./obj/main.o ./obj/functions.o
	$(CC) ./obj/main.o ./obj/functions.o -o ./bin/matrix $(CFLAGS)

./obj/main.o: ./src/main.c
	$(CC) -c ./src/main.c -o ./obj/main.o $(CFLAGS)

./obj/functions.o: ./src/functions.c
	$(CC) -c ./src/functions.c -o ./obj/functions.o $(CFLAGS)

clean:
	rm -f ./obj/*.o
	rm -f ./bin/matrix
