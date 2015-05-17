CC = mpicc
LD = mpicc
CFLAGS = -g
LFLAGS = -O0

ALL_SRC := $(wildcard *.c)
EXCLUDES = mpi_dyno.c pseudocode.c $(wildcard *test.c)
SRCS := $(filter-out $(EXCLUDES),$(ALL_SRC))
OBJS := $(patsubst %.c,%.o,$(SRCS))

.PHONY: all
all: dynogrid

.PHONY: tests
tests: list_test
	$<

dynogrid: $(OBJS)
	$(LD) $(LFLAGS) $^ -o $@ -lm

list_test: list_test.o list.o
	$(LD) $(LFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: run
run: dynogrid
	@mpirun -np 2 -hostfile hostfile ./dynogrid 1 2 1

.PHONY: debug
debug: dynogrid
	@idev -l h_rt=01:00:00

.PHONY: valgrind
valgrind: dynogrid
	@valgrind -v --leak-check=full ./dynogrid

clean:
	@/bin/rm -f *.d
	@/bin/rm -f idev*
	@/bin/rm -f *.o
	@/bin/rm -f dynogrid

-include $(OBJS:%.o=%.d)
