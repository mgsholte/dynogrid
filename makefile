CC = mpicc
CFLAGS = -g
LFLAGS = -O3

ALL_SRC := $(patsubst %.c,%.o,$(wildcard *.c))
EXCLUDES = 
OBJS := $(filter-out $(EXCLUDES),$(ALL_SRC))


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
	@rm -f batch_dynogrid
	@qsub jobscript
	@watch -n 10 qstat

.PHONY: valgrind
valgrind: dynogrid
	@valgrind -v --leak-check=full ./dynogrid

clean:
	@/bin/rm -f *.o
	@/bin/rm -f dynogrid
