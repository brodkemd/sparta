# Compiler and flags
CC = mpic++
ALLFLAGS = -DSPARTA_UNORDERED_MAP -DSPARTA_BIGBIG -g
ELMER_FLAGS = $(shell python3-config --cflags)
SPARTA_FLAGS = -O3
SHFLAGS = -fPIC
LDFLAGS = $(shell python3-config --ldflags --embed)
EXTRA_LDFLAGS = -O

EXE = sparta_fea
SOURCES := $(wildcard *.cpp) $(wildcard ../*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

# Implicit rule would not work as object files are in different folders
all: $(EXE)

configure:
	python3 ./shell_server_config/make_shell_server_config.py
	@echo ""
	python3 ./configure.py

$(EXE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $(EXE)

../%.o: ../%.cpp
	$(CC) $(ALLFLAGS) $(SPARTA_FLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(ALLFLAGS) $(ELMER_FLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm -f $(OBJECTS) $(EXE)
