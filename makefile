CC = mpicc
CFLAGS = -g
LFLAGS = -O3

ALL_SRC := $(patsubst %.c,%.o,$(wildcard *.c))
EXCLUDES = pseudocode.o balance.o
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

.PHONY: debug
debug: dynogrid
	@idev -l h_rt=02:00:00

.PHONY: valgrind
valgrind: dynogrid
	@valgrind -v --leak-check=full ./dynogrid

clean:
	@/bin/rm -f *.o
	@/bin/rm -f dynogrid
