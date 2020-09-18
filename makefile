FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o
	gcc main.o Divide_Into_Modularity_Groups.o IO.o Stack.o Divide_Into_Two.o Flag_Set.o Linked_List.o Leading_eigenpair.o modMat.o Spmat.o  -o cluster $(LIBS)
clean:
	rm -rf *.o cluster

main.o: main.c Divide_Into_Modularity_Groups.o IO.o Stack.o 
	gcc $(FLAGS) -c main.c

IO.o: IO.c IO.h
	gcc $(FLAGS) -c IO.c

Stack.o: Stack.c Stack.h
	gcc $(FLAGS) -c Stack.c	 

Divide_Into_Modularity_Groups.o: Divide_Into_Modularity_Groups.c Divide_Into_Modularity_Groups.h Divide_Into_Two.o modMat.o 
	gcc $(FLAGS) -c Divide_Into_Modularity_Groups.c

Divide_Into_Two.o: Divide_Into_Two.c Divide_Into_Two.h Leading_eigenpair.o Flag_Set.o Linked_List.o
	gcc $(FLAGS) -c Divide_Into_Two.c

Flag_Set.o: Flag_Set.c Flag_Set.h
	gcc $(FLAGS) -c Flag_Set.c
	
Linked_List.o: Linked_List.c Linked_List.h
	gcc $(FLAGS) -c Linked_List.c

Leading_eigenpair.o: Leading_eigenpair.c Leading_eigenpair.h modMat.o 
	gcc $(FLAGS) -c Leading_eigenpair.c

modMat.o: modMat.c modMat.h Spmat.o
	gcc $(FLAGS) -c modMat.c
	
Spmat.o: Spmat.c Spmat.h
	gcc $(FLAGS) -c Spmat.c

print_bin: print_bin.o
	gcc print_bin.o -o print_bin $(LIBS)

print_bin.o: print_bin.c
	gcc $(FLAGS) -c print_bin.c