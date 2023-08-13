CXX:=g++
CXXFLAGS:= -Ofast -Wall -std=c++20 -fopenmp

EXT = cpp

SOURCES := $(wildcard *.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),%.o,$(SOURCES))

all: test speed_test

debug: CXXFLAGS += -g
debug: test speed

test: main.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

speed_test:  speed_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	clear && rm test speed_test
