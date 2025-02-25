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

all: basic_analytics

basic_analytics: basic_analytics.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) basic_analytics.cpp -o basic_analytics

clean:
	rm basic_analytics
