FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o
	gcc main.o Divide_Into_Modularity_Groups.o IO.o Stack.o divide_into_two.o leading_eigenpair.o modMat.o spmat.o  -o cluster $(LIBS)
clean:
	rm -rf *.o run

main.o: main.c Divide_Into_Modularity_Groups.o IO.o Stack.o 
	gcc $(FLAGS) -c main.c

IO.o: IO.c IO.h
	gcc $(FLAGS) -c IO.c

Stack.o: Stack.c Stack.h
	gcc $(FLAGS) -c Stack.c	 

Divide_Into_Modularity_Groups.o: Divide_Into_Modularity_Groups.c Divide_Into_Modularity_Groups.h divide_into_two.o modMat.o 
	gcc $(FLAGS) -c Divide_Into_Modularity_Groups.c

divide_into_two.o: divide_into_two.c divide_into_two.h leading_eigenpair.o 
	gcc $(FLAGS) -c divide_into_two.c
	
leading_eigenpair.o: leading_eigenpair.c leading_eigenpair.h modMat.o 
	gcc $(FLAGS) -c leading_eigenpair.c

modMat.o: modMat.c modMat.h spmat.o
	gcc $(FLAGS) -c modMat.c
	
spmat.o: spmat.c spmat.h
	gcc $(FLAGS) -c spmat.c


