# Makefile

CXX = g++
CXXFLAGS = -fdiagnostics-color -std=c++17 -O3 -Wall -Werror -Wextra
INCLUDES = -I/usr/include/hdf5/serial
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5

# Targets
TARGETS = brute.out hnsw.out

# Build all targets
all: $(TARGETS)

# Compile brute_force_search.cpp
brute.out: brute_force_search.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@

# Compile hnsw.cpp
hnsw.out: hnsw.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@

# Clean up compiled files
clean:
	rm -f $(TARGETS)
