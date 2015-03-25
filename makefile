OBJS = main.o push.o
CFLAGS = -O3 -openmp -c
LFLAGS = -O3 -openmp
all:scan

CC=dyno

dyno: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o dynogrid

main.o: main.c
	$(CC) $(CFLAGS) main.c

push.o: push.c
	$(CC) $(CFLAGS) push.c

run:scan
	@./dynogrid

clean:
	@/bin/rm -rf *.o
