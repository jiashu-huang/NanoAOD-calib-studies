CXX := g++
CXXFLAGS := -O2 -Wall -Wextra -Wpedantic -std=c++17
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
GIT_COMMIT := $(shell git describe --always --dirty 2>/dev/null || echo unknown)

match-filter: MatchFilter.cc
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) -DGIT_COMMIT='"$(GIT_COMMIT)"' $< -o $@ $(ROOT_LIBS)

.PHONY: clean
clean:
	rm -f match-filter
