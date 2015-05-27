CC=gcc
CXX=g++

CFLAGS = -g -Wall -I./shared
#CFLAGS = -O3 -Wall -I./shared
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LIB = -lgsl -lgslcblas -lz -lpthread

SHARED_LIB = gen_func.cpp read_data.cpp threadpool.c bfgs.cpp



all: $(SHARED_LIB) parse_args EM ngsF-HMM
	$(CXX) $(DFLAGS) *.o $(LIB) -o ngsF-HMM



$(SHARED_LIB):
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared/$@

parse_args:
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

EM:
	$(CXX) $(CFLAGS) $(DFLAGS) -c EM.cpp

ngsF-HMM:
	$(CXX) $(CFLAGS) $(DFLAGS) -c ngsF-HMM.cpp

test:
	@cd examples/; bash ./test.sh 2> test.log; cd ../

clean:
	@rm -f *~ *.o ngsF-HMM examples/testF* examples/testF.log
