CXX:=g++
CXXFLAGS:= -O3 -Wall -std=c++20 -fopenmp

EXT = cpp

SOURCES := $(wildcard *.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),%.o,$(SOURCES))


all: test speed_test lin_alg

profile: CXXFLAGS:= -O3 -pg -Wall -std=c++20
profile: all

debug: CXXFLAGS:= -Og -fopenmp -g -Wall -std=c++20
debug: all

test: main.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

lin_alg: lin_alg.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

speed_test:  speed_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	clear && rm test speed_test lin_alg
