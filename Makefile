CXX:=g++
CXXFLAGS:= -Ofast -Wall -std=c++17 -fopenmp

EXT = cpp

SOURCES := $(wildcard *.$(EXT))
OBJECTS := $(patsubst $(SRC)/%.$(EXT),%.o,$(SOURCES))

all: run

debug: CXXFLAGS += -g
debug: run

run: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.$(EXT) Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	clear && rm run
