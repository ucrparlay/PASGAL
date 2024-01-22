ifdef GCC
CC = g++
else
CC = clang++
endif

CPPFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../../external/parlaylib/include/ -I../../external/ -I.. 

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

all: fast-bcc tarjan-vishkin hopcroft-tarjan

fast-bcc:	fast-bcc.cpp fast-bcc.h connectivity.h ldd.h spanning_forest.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) fast-bcc.cpp -o fast-bcc

tarjan-vishkin:	tarjan-vishkin.cpp tarjan-vishkin.h connectivity.h ldd.h spanning_forest.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) tarjan-vishkin.cpp -o tarjan-vishkin

hopcroft-tarjan:	hopcroft-tarjan.cpp hopcroft-tarjan.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) hopcroft-tarjan.cpp -o hopcroft-tarjan

clean:
	rm fast-bcc tarjan-vishkin hopcroft-tarjan
