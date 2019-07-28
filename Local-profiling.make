SHELL='bash'
#
# Bridges - PSC
#
# Intel Compilers are loaded by default
# You will need to specifically switch to GNU Modules
# With with `modules.sh`
#

CC = g++
MPCC = mpic++
OPENMP = -fopenmp
CFLAGS = -O3
PROFILING = -pg
LIBS =


TARGETS = serial pthreads openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o $(PROFILING)
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o $(PROFILING)
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -pthread pthreads.o common.o $(PROFILING)
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o $(PROFILING)
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o $(PROFILING)

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp $(PROFILING)
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp $(PROFILING)
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp $(PROFILING)
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp $(PROFILING)
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp $(PROFILING)
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp $(PROFILING)

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
