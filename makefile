FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main.o
	gcc main.o Divide_Into_Modularity_Groups.o IO.o Stack.o Divide_Into_Two.o FlagSet.o Linked_List.o leading_eigenpair.o modMat.o spmat.o  -o cluster $(LIBS)
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

Divide_Into_Two.o: Divide_Into_Two.c Divide_Into_Two.h leading_eigenpair.o FlagSet.o Linked_List.o
	gcc $(FLAGS) -c Divide_Into_Two.c

FlagSet.o: FlagSet.c FlagSet.h
	gcc $(FLAGS) -c FlagSet.c
	
Linked_List.o: Linked_List.c Linked_List.h
	gcc $(FLAGS) -c Linked_List.c

leading_eigenpair.o: leading_eigenpair.c leading_eigenpair.h modMat.o 
	gcc $(FLAGS) -c leading_eigenpair.c

modMat.o: modMat.c modMat.h spmat.o
	gcc $(FLAGS) -c modMat.c
	
spmat.o: spmat.c spmat.h
	gcc $(FLAGS) -c spmat.c

print_bin: print_bin.o
	gcc print_bin.o -o print_bin $(LIBS)

print_bin.o: print_bin.c
	gcc $(FLAGS) -c print_bin.c