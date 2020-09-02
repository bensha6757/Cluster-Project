#ifndef IO_MEM_ERRORS_H_
#define IO_MEM_ERRORS_H_

#include <stdio.h>
#include <stdlib.h>


#define VERIFY(expression, error)                                                           \
    if (!(expression)){                                                                      \
        printf("\n%s%s, %s", "Verification Failed: ", #expression, stringFromError(error));   \
        printf("\n");                                                                          \
        exit(error);                                                                               \
    }                                                                \

typedef enum error_t {
    SUCCESS,
    MISSING_ARG_ERROR,
    FILE_NOT_FOUND_ERROR,
	FILE_READ_ERROR,
	FILE_WRITE_ERROR,
	MEM_ALLOC_ERROR,
	INFINITE_LOOP_ERROR,
    OUT_OF_BOUNDS_ERROR,
    NULL_POINTER_ERROR
} ERROR;

static const char* stringFromError(enum error_t e)
{
    static const char *strings[] = { 
                                    "Program Completed Successfully",
                                    "Missing program argument error",
                                    "Error while opening file: not found or corrupted",
                                    "Error while reading the file", 
                                    "Error while writing to file", 
                                    "Memory allocation error", 
                                    "Infinite loop error", 
                                    "Out of bounds error", 
                                    "Null pointer error", 
                                    };
    return strings[e];
}

#endif /* IO_MEM_ERRORS_H_ */
