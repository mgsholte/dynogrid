OBJS = search.o bsearch.o
CFLAGS = -O3 -openmp -c
LFLAGS = -O3 -openmp
all:scan

CC=mpicxx

scan: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o dynogrid

bsearch.o: bsearch.cpp
	$(CC) $(CFLAGS) bsearch.cpp

search.o: search.cpp
	$(CC) $(CFLAGS) search.cpp

run:scan
	@./search 10

clean:
	@/bin/rm -rf *.o
