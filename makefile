CC = mpicc
CFLAGS = -g
LFLAGS = -O3

OBJS := $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all
all: dynogrid

dynogrid: $(OBJS)
	$(CC) $(LFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $^ -o $@

push.c: dynamics.h decs.h
	@touch $@

%.c: %.h decs.h
	@touch $@

.PHONY: run
run: dynogrid
	@./dynogrid

clean:
	@/bin/rm -rf *.o
