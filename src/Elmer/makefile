# Compiler and flags
CC = mpic++
CFLAGS = $(shell python3-config --cflags)
SHFLAGS = -fPIC
LDFLAGS = -O $(shell python3-config --ldflags --embed)

EXE = sparta_fea
SOURCES := $(wildcard *.cpp) $(wildcard ../*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

# Implicit rule would not work as object files are in different folders
all: $(EXE)

$(EXE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(EXE)
    # $(CC) $(LDFLAGS) $^ -o $@

../%.o: ../%.cpp
    # $(CC) $(CFLAGS) $(SHFLAGS) -I../ -c $< -o $@
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
    # $(CC) $(CFLAGS) $(SHFLAGS) -c $< -o $@
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(OBJECTS) $(EXE)