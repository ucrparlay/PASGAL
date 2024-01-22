ifdef GCC
CC = g++
else
CC = clang++
endif

CPPFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../../external/parlaylib/include/ -I..

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

all: scc tarjan

scc:	scc.cpp scc.h reach.h multi_reach.h resizable_table.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) scc.cpp -o scc

scc_seq: scc.cpp scc.h reach.h multi_reach.h resizable_table.h
	$(CC) -std=c++20 -Wall -Wextra -Werror -DPARLAY_SEQUENTIAL -O3 -mcx16 -march=native  $(INCLUDE_PATH) scc.cpp -o scc_seq
tarjan: tarjan.h tarjan.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) tarjan.cpp -o tarjan

clean:
	rm scc tarjan
