CC=gcc
CXX=g++

CFLAGS = -g -Wall
#CFLAGS = -O3 -Wall
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LIB = -lgsl -lgslcblas -lz -lpthread

all: ngsF-HMM


parse_args: parse_args.cpp ngsF-HMM.hpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

threadpool: threadpool.c threadpool.h
	$(CXX) $(CFLAGS) $(DFLAGS) -c threadpool.c

EM: EM.cpp ngsF-HMM.hpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c EM.cpp

shared: shared.cpp shared.hpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared.cpp

ngsF-HMM: ngsF-HMM.cpp parse_args threadpool EM shared
	$(CXX) $(CFLAGS) $(DFLAGS) ngsF-HMM.cpp parse_args.o threadpool.o EM.o shared.o $(LIB) -o ngsF-HMM

test:
	@cd examples/; bash ./test.sh 2> test.log; cd ../

clean:
	@rm -f *~ *.o ngsF-HMM examples/testF*
