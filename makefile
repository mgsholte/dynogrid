CC = mpicc
CFLAGS = -O3 -openmp -c
LFLAGS = -O3 -openmp

OBJS = main.o push.o

.PHONY: all
all: dynogrid

dynogrid: $(OBJS)
	$(CC) $(LFLAGS) $^ -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $^

.PHONY: run
run: dynogrid
	@./dynogrid

clean:
	@/bin/rm -rf *.o
