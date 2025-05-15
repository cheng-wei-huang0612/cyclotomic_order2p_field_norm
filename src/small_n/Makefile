CC = gcc
CFLAGS = -O2 -Wall -std=c99
LDLIBS = -lgmp

all: norm_schoolbook norm_ntt

norm_schoolbook: norm_schoolbook.c utils.c utils.h
	$(CC) $(CFLAGS) norm_schoolbook.c utils.c -o norm_schoolbook $(LDLIBS)

norm_ntt: norm_ntt.c utils.c utils.h
	$(CC) $(CFLAGS) norm_ntt.c utils.c -o norm_ntt $(LDLIBS)

clean:
	rm -f norm_schoolbook norm_ntt

