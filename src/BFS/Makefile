ifdef GCC
CC = g++
else
CC = clang++
endif

CPPFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../../external/parlaylib/include/ -I../ -I../../external/

ifdef CILKPLUS
CC = clang++
CPPFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CPPFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CPPFLAGS += -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -pthread
endif

ifdef DEBUG
CPPFLAGS += -DDEBUG -O0 -g -fsanitize=address
else ifdef PERF
CC = g++
CPPFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CPPFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL -g
else
CPPFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CPPFLAGS += -DPARLAY_USE_STD_ALLOC
endif

ifdef STAT
CPPFLAGS += -DSTAT
endif

all: bfs seq-bfs bfs_test bfs_delta bfs_delta_winning_tree

bfs:	bfs.cpp bfs.h seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs.cpp -o bfs

bfs_chunk:	bfs.cpp bfs.h seq-bfs.h ../chunkbag.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) -DBAG_CHUNK bfs.cpp -o bfs_chunk

bfs_delta_bench:	bfs_delta_bench.cpp bfs_delta.h ../chunkbag.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs_delta_bench.cpp -o bfs_delta_bench

# Sanitizer builds (separate binaries; asserts stay enabled via -UNDEBUG)
SAN_FLAGS = -std=c++20 -Wall -Wextra -O1 -g -pthread $(INCLUDE_PATH)

sanitizers: bfs_tsan bfs_chunk_tsan bfs_asan bfs_chunk_asan

bfs_tsan:	bfs.cpp bfs.h seq-bfs.h
	clang++ $(SAN_FLAGS) -fsanitize=thread bfs.cpp -o bfs_tsan

bfs_chunk_tsan:	bfs.cpp bfs.h seq-bfs.h ../chunkbag.h
	clang++ $(SAN_FLAGS) -fsanitize=thread -DBAG_CHUNK bfs.cpp -o bfs_chunk_tsan

bfs_asan:	bfs.cpp bfs.h seq-bfs.h
	clang++ $(SAN_FLAGS) -fsanitize=address bfs.cpp -o bfs_asan

bfs_chunk_asan:	bfs.cpp bfs.h seq-bfs.h ../chunkbag.h
	clang++ $(SAN_FLAGS) -fsanitize=address -DBAG_CHUNK bfs.cpp -o bfs_chunk_asan

seq-bfs:	seq-bfs.cpp seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) seq-bfs.cpp -o seq-bfs

bfs_test:	bfs_test.cpp bfs.h seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs_test.cpp -o bfs_test

bfs_delta:	bfs_delta.cpp bfs_delta.h seq-bfs.h ../chunkbag.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs_delta.cpp -o bfs_delta

bfs_delta_winning_tree:	bfs_delta_winning_tree.cpp bfs_delta_winning_tree.h seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs_delta_winning_tree.cpp -o bfs_delta_winning_tree

reachability:	reachability.cpp reachability.h seq-bfs.h ../chunkbag.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) reachability.cpp -o reachability

reachability_bench:	reachability_bench.cpp reachability.h ../chunkbag.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) reachability_bench.cpp -o reachability_bench

clean:
	rm bfs seq-bfs bfs_test
