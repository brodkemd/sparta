CFLAGS = -O
CC = g++
SRC = toml.cpp main.cpp fixfea.cpp
OBJ = $(SRC:.cpp = .o)

NumTest: $(OBJ)
	$(CC) $(CFLAGS) -o run $(OBJ)

clean:
	rm run *.o