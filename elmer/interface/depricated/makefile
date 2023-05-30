run: build/make_sif.o build/elmer.o
	g++ build/make_sif.o build/elmer.o -o run

build/make_sif.o:
	g++ -c make_sif.cpp -o build/make_sif.o

build/elmer.o:
	g++ -c elmer.cpp -o build/elmer.o

# util.o:
# 	g++ -c tools/util.cpp -o build/util.o

clean:
	rm build/*.o run