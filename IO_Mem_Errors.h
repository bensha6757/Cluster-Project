#ifndef IO_Mem_Errors_H_
#define IO_Mem_Errors_H_

#include <stdio.h>
#include <stdlib.h>

/******************************************************************
 * a module in charge of exiting the program in case of an error, *
 * and supplying a message explaining the reason                  *
 * ****************************************************************/ 

typedef enum exit_code_t {
    SUCCESS,
    MISSING_ARG_ERROR,
    FILE_NOT_FOUND_ERROR,
	FILE_READ_ERROR,
	FILE_WRITE_ERROR,
	MEM_ALLOC_ERROR,
	INFINITE_LOOP_ERROR,
    OUT_OF_BOUNDS_ERROR,
    NULL_POINTER_ERROR,
    DIVISION_BY_ZERO
} exit_code;

/* assertion check */
#define VERIFY(expression, error)                                                           \
    if (!(expression)){                                                                      \
        printf("\n%s%s, %s", "Verification Failed: ", #expression, stringFromError(error));   \
        printf("\n");                                                                          \
        exit(error);                                                                            \
    }                                                                                            \

/* print the relevent error message */
static const char* stringFromError(exit_code e)
{
    static const char *strings[] = 
    { 
        "Program Completed Successfully",
        "Missing program argument error",
        "Error while opening file: not found or corrupted",
        "Error while reading the file", 
        "Error while writing to file", 
        "Memory allocation error", 
        "Infinite loop error", 
        "Out of bounds error", 
        "Null pointer error",
        "Division by zero error" 
    };
    return strings[e];
}

#endif /* IO_Mem_Errors_H_ */
