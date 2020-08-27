#ifndef IO_MEM_ERRORS_H_
#define IO_MEM_ERRORS_H_

#include <stdio.h>
#include <stdlib.h>


#define VERIFY(expression, error)                                                           \
    if (!(expression)){                                                                      \
        printf("\n%s%s, %s", "Verification Failed: ", #expression, stringFromError(error));   \
        exit(1);                                                                               \
    }                                                                \

typedef enum error_t {
	FILE_READ_ERROR,
	FILE_WRITE_ERROR,
	MEM_ALLOC_ERROR,
	INFINITE_LOOP_ERROR,
    OUT_OF_BOUNDS_ERROR,
    NULL_POINTER_ERROR
} ERROR;

static const char* stringFromError(enum error_t e)
{
    static const char *strings[] = { "Error while reading the file", "Error while writing to file", "Memory allocation error", "Infinite loop error", "Out of bounds error", "Null pointer error"};
    return strings[e];
}

#endif /* IO_MEM_ERRORS_H_ */
