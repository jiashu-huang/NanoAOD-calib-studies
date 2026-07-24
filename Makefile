CXX := g++
CXXFLAGS := -O2 -Wall -Wextra -Wpedantic -std=c++17
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

calib-study: CalibStudy.cc
	$(CXX) $(CXXFLAGS) $(ROOT_CFLAGS) $< -o $@ $(ROOT_LIBS)

.PHONY: clean
clean:
	rm -f calib-study
