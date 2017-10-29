
exec: FORCE main.o
	gcc -fopenmp -o exec main.o
	./exec

main.o: main.c
	gcc -c main.c -fopenmp -o main.o

FORCE:
