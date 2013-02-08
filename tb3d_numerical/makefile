#CC=mpiCC
CC=mpiCC -O2 -Wno-unused-result
tb.out: grid.o multigrid.o mesh.o io.o main.o
	$(CC) main.o grid.o multigrid.o mesh.o io.o -o tb.out -lm -lboost_mpi -lboost_serialization -lboost_program_options -lmpreal -lmpfr -lgmp

grid.o: grid.cpp grid.h
	$(CC) -c grid.cpp 
multigrid.o: multigrid.cpp multigrid.h
	$(CC) -c multigrid.cpp 
mesh.o: mesh.cpp mesh.h
	$(CC) -c mesh.cpp 
io.o: io.cpp io.h
	$(CC) -c io.cpp 
main.o: main.cpp multigrid.o mesh.o grid.o
	$(CC) -c main.cpp 
clean:
	rm *.o *.out

