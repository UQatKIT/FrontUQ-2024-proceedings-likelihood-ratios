# Compiler
CC = g++
# Compiler options
CFLAGS = -g -Wall -O3 -DNDEBUG -I /usr/local/include/eigen3 -fopenmp -I /home/jj6547/.local/eigen3 -fopenmp


# Target executable
TARGET = build/likelihoodRatioMC

# all .cpp files in this directory are my sources
SOURCES = $(wildcard *.cpp)

# .o files depend upon the .cpp files
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

build/likelihoodRatioMC: likelihoodRatioMC.o heat_equation_mc.o heat_equation_grid.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) $(LIB) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)
