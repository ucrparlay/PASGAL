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
CPPFLAGS += -DDEBUG -Og -g
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

all: bfs seq-bfs bfs_test

bfs:	bfs.cpp bfs.h seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs.cpp -o bfs

seq-bfs:	seq-bfs.cpp seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) seq-bfs.cpp -o seq-bfs

bfs_test:	bfs_test.cpp bfs.h seq-bfs.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) bfs_test.cpp -o bfs_test

clean:
	rm bfs seq-bfs bfs_test
