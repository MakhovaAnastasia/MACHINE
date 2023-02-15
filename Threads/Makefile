all: ref

		
ref: main.o solution.o input.o discrepancy.o
	g++ main.o solution.o input.o discrepancy.o -o ref -pthread -O3

main.o: main.cpp
	g++ -c main.cpp -pthread -O3

solution.o: solution.cpp
	g++ -c solution.cpp -pthread -O3

input.o: input.cpp
	g++ -c input.cpp -pthread -O3

discrepancy.o: discrepancy.cpp
	g++ -c discrepancy.cpp -pthread -O3

clean:
		rm -rf *.o ref
		
cleano:
		rm -rf *.o
