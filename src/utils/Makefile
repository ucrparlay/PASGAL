ifdef GCC
CC = g++
else
CC = clang++
endif

CPPFLAGS = -std=c++20 -Wall -Wextra -Werror

INCLUDE_PATH = -I../external/parlaylib/include/ -I.. -I../../external/

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
CPPFLAGS += -DDEBUG -Og
else ifdef PERF
CC = g++
CPPFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CPPFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CPPFLAGS += -DPARLAY_USE_STD_ALLOC
endif

all: symmetrize pbbs2bin generate_random_graph generate_grid_graph get_diameter

symmetrize:	symmetrize.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) symmetrize.cpp -o symmetrize

pbbs2bin:	pbbs2bin.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) pbbs2bin.cpp -o pbbs2bin

validator:	validator.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) validator.cpp -o validator

generate_random_graph:	generate_random_graph.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) generate_random_graph.cpp -o generate_random_graph

generate_grid_graph:	generate_grid_graph.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) generate_grid_graph.cpp -o generate_grid_graph

get_diameter:	get_diameter.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) get_diameter.cpp -o get_diameter

clean:
	rm symmetrize pbbs2bin validator generate_random_graph generate_grid_graph get_diameter
