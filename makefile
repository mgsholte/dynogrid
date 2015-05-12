CC = mpicc
LD = mpicc
CFLAGS = -g -MMD -Wcheck -Werror
LFLAGS = -O3

ALL_SRC := $(wildcard *.c)
EXCLUDES = pseudocode balance
SRCS := $(filter-out $(addsuffix .c,$(EXCLUDES)),$(ALL_SRC))
OBJS := $(patsubst %.c,%.o,$(SRCS))

.PHONY: all
all: dynogrid

dynogrid: $(OBJS)
	$(LD) $(LFLAGS) $^ -o $@ -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY: run
run: dynogrid
	@rm -f batch_dynogrid
	@qsub jobscript
	@watch -n 10 qstat

.PHONY: debug
debug: dynogrid
	@idev -l h_rt=01:00:00

.PHONY: valgrind
valgrind: dynogrid
	@valgrind -v --leak-check=full ./dynogrid

clean:
	@/bin/rm -f *.d
	@/bin/rm -f *.o
	@/bin/rm -f dynogrid

-include $(OBJS:%.o=%.d)
