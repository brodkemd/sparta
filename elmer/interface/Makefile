run: build/elmer_to_sparta_temperature.o build/elmer_grid_to_sparta_surf.o build/main.o build/util.o
	g++ build/elmer_to_sparta_temperature.o build/elmer_grid_to_sparta_surf.o build/main.o build/util.o -o run

build/elmer_to_sparta_temperature.o:
	g++ -c elmer_to_sparta_temperature.cpp -o build/elmer_to_sparta_temperature.o

build/elmer_grid_to_sparta_surf.o:
	g++ -c elmer_grid_to_sparta_surf.cpp -o build/elmer_grid_to_sparta_surf.o

build/main.o:
	g++ -c main.cpp -o build/main.o

build/util.o:
	g++ -c tools/util.cpp -o build/util.o

clean:
	rm build/*.o run