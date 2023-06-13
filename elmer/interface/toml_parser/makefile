run: out.o main.o handles.o
	g++ -o run out.o main.o handles.o

main.o: main.cpp
	g++ -o main.o -c main.cpp

out.o: out.cpp
	g++ -o out.o -c out.cpp 

handles.o: handles.cpp
	g++ -o handles.o -c handles.cpp 


clean:
	rm *.o run