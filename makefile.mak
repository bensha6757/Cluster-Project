FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o Divide_Into_Modularity_Groups.o IO.o Stack.o Types.h
	gcc main.o -o cluster $(LIBS)
clean:
	rm -rf *.o run

main.o: main.c io_mem_errors.h
	gcc $(FLAGS) -c main.c

IO.o: IO.c IO.h io_mem_errors.h
		gcc $(FLAGS) -C  Stack.c

Stack.o: Stack.c Stack.h io_mem_errors.h
	gcc $(FLAGS) -C  Stack.c	 

Divide_Into_Modularity_Groups.o: Divide_Into_Modularity_Groups.c Divide_Into_Modularity_Groups.h divide_into_two.o io_mem_errors.h
	gcc $(FLAGS) -C  Divide_Into_Modularity_Groups.c

divide_into_two.o: divide_into_two.c divide_into_two.h leading_eigenpair.o io_mem_errors.h
	gcc $(FLAGS) -C divide_into_two.c
	
leading_eigenpair.o: leading_eigenpair.c leading_eigenpair.h modMat.o io_mem_errors.h
	gcc $(FLAGS) -C leading_eigenpair.c

modMat.o: modMat.c modMat.h spmat.o io_mem_errors.h
	gcc $(FLAGS) -C modMat.c
	
spmat.o: spmat.c spmat.h io_mem_errors.h
	gcc $(FLAGS) -c spmat.c


