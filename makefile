CC = mpicc
CFLAGS = -O3 -openmp
LFLAGS = -O3 -openmp

OBJS := $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all
all: dynogrid

dynogrid: $(OBJS)
	$(CC) $(LFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $^ -o $@

push.c: dynamics.h decs.h

%.c: %.h decs.h

.PHONY: run
run: dynogrid
	@./dynogrid

clean:
	@/bin/rm -rf *.o
