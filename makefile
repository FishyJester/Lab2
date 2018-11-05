CC = gcc
CFLAGS = -lm -lpthread

.PHONY : clean

newton : newton.c
	$(CC) -o $@ $(CFLAGS) newton.c

clean :
	rm -f newton
