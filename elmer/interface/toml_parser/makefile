CFLAGS = -O
CC = g++
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp = .o)

NumTest: $(OBJ)
	$(CC) $(CFLAGS) -o run $(OBJ)

clean:
	rm run *.o