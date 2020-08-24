/*
 * io_mem_errors.h
 *
 *  Created on: 22 ����� 2020
 *      Author: ��
 */

#ifndef IO_MEM_ERRORS_H_
#define IO_MEM_ERRORS_H_

typedef enum error_t {
	FILE_READ_ERROR,
	FILE_WRITE_ERROR,
	MEM_ALLOC_ERROR,
	INFINITE_LOOP_ERROR,
    OUT_OF_BOUNDS_ERROR,
    NULL_POINTER_ERROR
} ERROR;

#define verify(expression, msg)                                         \
    if (!expression){                                                   \
        printf("\n%s%s", "Verification Failed: ", #expression);         \
        exit(msg);                                                      \
    }                                                                   \



#endif /* IO_MEM_ERRORS_H_ */
