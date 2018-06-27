
CFLAGS=-O3 -Wall -Wextra -Werror -msse4.1 -mssse3 -msha

all: sha256qc

clean:
	rm -f sha256i.o

sha256i: sha256i.o
	$(CC) $(CFLAGS) $< -o $@

sha256qc: sha256i.o sha256.hs Makefile
	ghc -O --make sha256.hs sha256i.o -o $@
