CXX:=g++
CXXFLAGS:= -O3 -Wall -std=c++20 -fopenmp -funroll-loops

EXT = cpp

SOURCES := $(wildcard *.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),%.o,$(SOURCES))


all: test speed_test

profile: CXXFLAGS:= -O3 -Wall -std=c++20
profile: speed_test

debug: CXXFLAGS:= -Og -g -Wall -std=c++20
debug: all

single: CXXFLAGS:= -O3 -Wall -std=c++20 -funroll-loops
single: all

test: tests/test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

speed_test:  tests/speed_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm test speed_test
	clear
