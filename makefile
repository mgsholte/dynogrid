CC = clang
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

.PHONY: valgrind
valgrind: dynogrid
	@valgrind -v --leak-check=full ./dynogrid

clean:
	@/bin/rm -f *.o
	@/bin/rm -f dynogrid
