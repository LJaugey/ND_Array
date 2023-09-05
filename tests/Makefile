CXX:=g++
CXXFLAGS:= -O3 -Wall -std=c++20 -fopenmp

EXT = cpp

SOURCES := $(wildcard *.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),%.o,$(SOURCES))

all: run

debug: CXXFLAGS:= -O3 -g -Wall -std=c++20
debug: all

test: test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: test
	./test

clean:
	clear && rm test